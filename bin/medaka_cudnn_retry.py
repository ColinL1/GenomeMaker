#!/usr/bin/env python3
"""`medaka` with a per-batch cuDNN fallback.

usage: medaka_cudnn_retry.py <medaka subcommand and arguments>

medaka 2.2.1 runs short contigs (< chunk_len) one at a time, unchunked. A contig whose
pileup has ~115k columns (insertions included) is rejected by the cuDNN GRU
(CUDNN_STATUS_NOT_SUPPORTED) and medaka inference exits. That batch is rerun with
cuDNN disabled (torch's native GRU kernel); all other batches still run through cuDNN.
On CPU the fallback never triggers. Retries are counted on stderr.
"""
import sys

import medaka.models
import torch
from medaka.medaka import main

_predict_on_batch = medaka.models.TorchModel.predict_on_batch
n_retry = 0


def predict_on_batch(self, batch):
    global n_retry
    try:
        return _predict_on_batch(self, batch)
    except RuntimeError as e:
        if 'CUDNN_STATUS_NOT_SUPPORTED' not in str(e):
            raise
        n_retry += 1
        shape = tuple(self.get_model_input_features(batch).shape)
        print(f"cuDNN retry {n_retry}: input shape {shape}", file=sys.stderr, flush=True)
        with torch.backends.cudnn.flags(enabled=False):
            return _predict_on_batch(self, batch)


medaka.models.TorchModel.predict_on_batch = predict_on_batch
sys.argv = ['medaka'] + sys.argv[1:]
rc = main()
print(f"cuDNN retries: {n_retry}", file=sys.stderr)
sys.exit(rc)
