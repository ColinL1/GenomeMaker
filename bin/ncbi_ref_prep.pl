#!/usr/bin/env perl
# Prepare an NCBI genome + annotation for `funannotate annotate`
# (PREP_REF_ANNOTATION in modules/10_comparison.nf).
#
#   ncbi_ref_prep.pl fasta <genomic.fna> <out.fa> <out.seqmap.tsv>
#       Trims headers to the accession and shortens names longer than 16
#       characters by dropping the ".N" version suffix (funannotate/tbl2asn
#       limit). Organelle sequences are dropped: they use another genetic code
#       and the nuclear-only Porites_sp annotation has no counterpart.
#       seqmap columns: original, new name or DROPPED, length, reason.
#
#   ncbi_ref_prep.pl gff <in.gff3> <seqmap.tsv>
#       Reduces a GFF3 to the bare gene models funannotate should re-annotate:
#       gene -> mRNA|tRNA|rRNA -> exon/CDS. Dropped: pseudogenes, mRNAs without
#       CDS, every other feature type, features on dropped sequences. Attributes
#       are cut down to ID/Parent (product kept for tRNA/rRNA) so that NCBI gene
#       names, products and Dbxrefs do not leak into the new annotation, and IDs
#       are restricted to [A-Za-z0-9._-]. funannotate uses the mRNA ID verbatim
#       as transcript ID, which is what keeps gffread proteins in step with it.
#
# The gff mode writes GFF3 to stdout and a summary to stderr.
use strict;
use warnings;

my $mode = shift @ARGV or die "usage: $0 fasta <fna> <out.fa> <map.tsv> | gff <gff3> <map.tsv>\n";

if ($mode eq 'fasta') {
    my ($fna, $out, $map) = @ARGV;
    die "usage: $0 fasta <fna> <out.fa> <map.tsv>\n" unless defined $map;
    open my $in, '<', $fna or die "$fna: $!\n";
    open my $fa, '>', $out or die "$out: $!\n";
    open my $mp, '>', $map or die "$map: $!\n";
    my (%seen, $orig, $new, $reason, $len);
    my $flush = sub {
        print $mp join("\t", $orig, $new, $len, $reason), "\n" if defined $orig;
    };
    while (my $line = <$in>) {
        if ($line =~ /^>(\S+)\s*(.*)/) {
            $flush->();
            ($orig, my $desc) = ($1, $2);
            ($new, $reason, $len) = ($orig, 'kept', 0);
            if ($desc =~ /mitochondri|chloroplast|plastid/i) {
                ($new, $reason) = ('DROPPED', 'organelle');
                next;
            }
            if (length($new) > 16) {
                $new =~ s/\.\d+$//;
                $reason = 'version suffix removed';
            }
            die "sequence name still longer than 16 characters: $orig\n" if length($new) > 16;
            die "duplicate sequence name after renaming: $new\n" if $seen{$new}++;
            print $fa ">$new\n";
        } else {
            next unless defined $new;
            (my $seq = $line) =~ s/\s+//g;
            $len += length $seq;
            print $fa $line unless $new eq 'DROPPED';
        }
    }
    $flush->();
}
elsif ($mode eq 'gff') {
    my ($gff, $map) = @ARGV;
    die "usage: $0 gff <gff3> <map.tsv>\n" unless defined $map;

    my %seqname;
    open my $mp, '<', $map or die "$map: $!\n";
    while (<$mp>) {
        chomp;
        my ($orig, $new) = split /\t/;
        $seqname{$orig} = $new;
    }

    my %tx_type = map { $_ => 1 } qw(mRNA tRNA rRNA);
    my (@rows, %type, %has_cds);
    open my $in, '<', $gff or die "$gff: $!\n";
    while (my $line = <$in>) {
        next if $line =~ /^#/ || $line !~ /\t/;
        chomp $line;
        my @f = split /\t/, $line, 9;
        next unless $f[2] =~ /^(gene|mRNA|tRNA|rRNA|exon|CDS)$/;
        my $seq = $seqname{$f[0]};
        die "sequence $f[0] is in the GFF but not in the FASTA\n" unless defined $seq;
        next if $seq eq 'DROPPED';
        $f[0] = $seq;
        my ($id)      = $f[8] =~ /(?:^|;)ID=([^;]+)/;
        my ($parent)  = $f[8] =~ /(?:^|;)Parent=([^;]+)/;
        my ($product) = $f[8] =~ /(?:^|;)product=([^;]+)/;
        my @parents = defined $parent ? split(/,/, $parent) : ();
        $type{$id} = $f[2] if defined $id && $f[2] !~ /^(exon|CDS)$/;
        if ($f[2] eq 'CDS') { $has_cds{$_} = 1 for @parents }
        push @rows, [\@f, $id, \@parents, $product];
    }

    # a transcript is kept when its parent is a real gene (not a pseudogene) and,
    # for mRNA, it has CDS; tRNA/rRNA need the product funannotate writes out
    my (%keep_tx, %keep_gene, %n);
    for my $r (@rows) {
        my ($f, $id, $parents, $product) = @$r;
        next unless $tx_type{$f->[2]} && defined $id && @$parents == 1;
        my $gene = $parents->[0];
        if (($type{$gene} // '') ne 'gene')               { $n{"$f->[2] not under a gene"}++; next }
        if ($f->[2] eq 'mRNA' && !$has_cds{$id})          { $n{'mRNA without CDS'}++;       next }
        if ($f->[2] ne 'mRNA' && !defined $product)       { $n{"$f->[2] without product"}++; next }
        $keep_tx{$id}     = 1;
        $keep_gene{$gene} = 1;
    }

    my $clean = sub { (my $x = shift) =~ s/[^A-Za-z0-9._-]/_/g; $x };
    my %clean_seen;
    for my $id (keys %keep_tx, keys %keep_gene) {
        my $c = $clean->($id);
        die "IDs collide after sanitising: $id and $clean_seen{$c}\n"
            if exists $clean_seen{$c} && $clean_seen{$c} ne $id;
        $clean_seen{$c} = $id;
    }

    print "##gff-version 3\n";
    my $sub = 0;
    for my $r (@rows) {
        my ($f, $id, $parents, $product) = @$r;
        my $t = $f->[2];
        if ($t eq 'gene') {
            next unless defined $id && $keep_gene{$id};
            print join("\t", @$f[0..7], 'ID=' . $clean->($id)), "\n";
            $n{gene}++;
        } elsif ($tx_type{$t}) {
            next unless defined $id && $keep_tx{$id};
            my $attr = 'ID=' . $clean->($id) . ';Parent=' . $clean->($parents->[0]);
            $attr .= ";product=$product" if $t ne 'mRNA';
            print join("\t", @$f[0..7], $attr), "\n";
            $n{$t}++;
        } else {
            for my $p (grep { $keep_tx{$_} } @$parents) {
                $sub++;
                print join("\t", @$f[0..7], "ID=" . lc($t) . "_$sub;Parent=" . $clean->($p)), "\n";
            }
        }
    }
    print STDERR "$_\t$n{$_}\n" for sort keys %n;
    die "no protein-coding gene models left after filtering\n" unless $n{mRNA};
}
else {
    die "unknown mode '$mode'\n";
}
