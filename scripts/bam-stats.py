#!/usr/bin/env python

import sys
import numpy as np
import pandas as pd
import pysam
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns
from functools import partial

tumor = pysam.AlignmentFile(snakemake.input[0], "rb")
normal = pysam.AlignmentFile(snakemake.input[1], "rb")

softclips = []

for i, rec in enumerate(normal):
    if rec.is_supplementary or rec.is_unmapped:
        continue
    is_first_read = rec.pos < rec.mpos
    get_clip = lambda c: c[1] if c[0] == 4 else None
    clip_left = get_clip(rec.cigartuples[0])
    if clip_left is not None:
        softclips.append([clip_left, True, is_first_read])
    clip_right = get_clip(rec.cigartuples[-1])
    if clip_right is not None:
        softclips.append([clip_right, False, is_first_read])
    if i == 10000000:
        break

softclips = pd.DataFrame(softclips, columns=["len", "left", "first_in_pair"])

def plot(*args, **kwargs):
    softclips = args[0]
    plt.hist(softclips, normed=True)
    q95 = np.percentile(softclips, 99)
    plt.plot([q95, q95], [0, 1.0], "--k")
    m = max(softclips)
    plt.plot([m, m], [0, 1.0], ":k")
    plt.text(m, 1, "max={}".format(m), horizontalalignment="right", verticalalignment="top")


g = sns.FacetGrid(softclips, col="left", row="first_in_pair")
g = g.map(plot, "len")

plt.savefig(snakemake.output[0])
