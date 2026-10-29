#!/usr/bin/env perl
# Prepare non-coding RNA annotations for the AGAT merge in MERGE_ANNOTATIONS.
#
#   ncrna_to_gff3.pl trnascan <tRNAscan-SE.gff>
#       Adds one explicit `gene` per tRNA. tRNAscan-SE's GFF has no gene
#       features, so AGAT invents genes and groups some tRNAs under one gene;
#       funannotate then keeps one tRNA per gene (302 tRNAs were lost).
#       tRNA pseudogenes pass through unchanged: `pseudogene` is already a
#       top-level feature for AGAT, which drops a gene put above it and leaves
#       a dangling Parent (fails gt gff3validator).
#
#   ncrna_to_gff3.pl rfam <cmscan.tblout> <Rfam family.txt>
#       Converts Infernal cmscan --fmt 2 hits to gene/ncRNA|rRNA/exon GFF3.
#       Kept: significant (inclusion '!'), non-overlapping ('*' or '^' in olp)
#       hits of Rfam families typed "Gene" in family.txt, except tRNA (tRNAscan-SE
#       covers tRNAs) and ribozymes. Cis-regulatory and other non-gene types are
#       dropped. rRNA families become `rRNA`, the rest `ncRNA` with ncRNA_class.
#
# Writes GFF3 to stdout.
use strict;
use warnings;

my $mode = shift @ARGV or die "usage: $0 trnascan <gff> | rfam <tblout> <family.txt>\n";

if ($mode eq 'trnascan') {
    my $gff = shift @ARGV or die "missing tRNAscan GFF\n";
    open my $in, '<', $gff or die "$gff: $!\n";
    print "##gff-version 3\n";
    while (my $line = <$in>) {
        next if $line =~ /^##gff-version/;
        if ($line =~ /^#/ || $line !~ /\t/) { print $line; next; }
        chomp $line;
        my @f = split /\t/, $line, 9;
        if ($f[2] eq 'tRNA') {
            my ($id) = $f[8] =~ /ID=([^;]+)/;
            print join("\t", @f[0..1], 'gene', @f[3..7], "ID=$id.gene"), "\n";
            (my $attr = $f[8]) =~ s/;?$//;
            print join("\t", @f[0..7], "$attr;Parent=$id.gene"), "\n";
        } else {
            print "$line\n";
        }
    }
    exit 0;
}

if ($mode eq 'rfam') {
    my ($tbl, $famfile) = @ARGV;
    die "missing tblout or family.txt\n" unless $tbl && $famfile;
    # family.txt (Rfam database dump): col1 accession, col2 id, col4 description,
    # col19 type (e.g. "Gene; snRNA; snoRNA; CD-box;", "Cis-reg;").
    my (%type, %desc);
    open my $ff, '<', $famfile or die "$famfile: $!\n";
    while (<$ff>) {
        chomp;
        my @c = split /\t/;
        next unless @c >= 19;
        $type{$c[0]} = $c[18];
        $desc{$c[0]} = $c[3];
    }
    open my $in, '<', $tbl or die "$tbl: $!\n";
    print "##gff-version 3\n";
    my $n = 0;
    while (<$in>) {
        next if /^#/;
        my @c = split ' ';
        next unless @c >= 20;
        # --fmt 2 columns (1-based): 2 target name, 3 accession, 4 query (seqid),
        # 10 seq from, 11 seq to, 12 strand, 17 score, 18 E-value, 19 inc, 20 olp
        my ($name, $acc, $seqid, $from, $to, $strand, $score, $evalue, $inc, $olp) =
            @c[1, 2, 3, 9, 10, 11, 16, 17, 18, 19];
        next unless $inc eq '!' && ($olp eq '*' || $olp eq '^');
        my $t = $type{$acc} // '';
        next unless $t =~ /^Gene/;
        next if $t =~ /tRNA;/ || $t =~ /ribozyme/;
        my ($start, $end) = $from <= $to ? ($from, $to) : ($to, $from);
        my $feat  = $t =~ /rRNA;/ ? 'rRNA' : 'ncRNA';
        my $class = $t =~ /snoRNA/            ? 'snoRNA'
                  : $t =~ /snRNA/             ? 'snRNA'
                  : $t =~ /miRNA/             ? 'miRNA'
                  : $name =~ /SRP/i           ? 'SRP_RNA'
                  : $name =~ /RNaseP/i        ? 'RNase_P_RNA'
                  : $name =~ /RNase_MRP/i     ? 'RNase_MRP_RNA'
                  : $name =~ /Vault/i         ? 'vault_RNA'
                  : $name =~ /^Y_RNA/i        ? 'Y_RNA'
                  : $name =~ /^7SK/i          ? 'other'
                  :                             'other';
        $n++;
        my $id   = sprintf "%s.rfam%d", $seqid, $n;
        (my $product = $desc{$acc} // $name) =~ s/[;=,]/ /g;
        my $attr = "ID=$id;Parent=$id.gene;Name=$name;product=$product;"
                 . "Dbxref=RFAM:$acc;rfam_score=$score;rfam_evalue=$evalue";
        $attr .= ";ncRNA_class=$class" if $feat eq 'ncRNA';
        print join("\t", $seqid, 'Infernal', 'gene', $start, $end, $score, $strand, '.', "ID=$id.gene;Name=$name"), "\n";
        print join("\t", $seqid, 'Infernal', $feat, $start, $end, $score, $strand, '.', $attr), "\n";
        print join("\t", $seqid, 'Infernal', 'exon', $start, $end, $score, $strand, '.', "ID=$id.exon1;Parent=$id"), "\n";
    }
    exit 0;
}

die "unknown mode '$mode'\n";
