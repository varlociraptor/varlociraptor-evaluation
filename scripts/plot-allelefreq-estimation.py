from itertools import product
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import common
import numpy as np


MIN_CALLS = 10


minlen = int(snakemake.wildcards.minlen)
maxlen = int(snakemake.wildcards.maxlen)
vartype = snakemake.wildcards.vartype

truth = common.load_variants(snakemake.input.truth, minlen, maxlen, vartype=vartype)
colors = common.get_colors(snakemake.config)

def plot(calls, colors):
    calls = calls[calls.is_tp]
    true_af = truth.loc[calls.MATCHING].reset_index().TAF
    calls = calls.reset_index()
    calls["error"] = calls.CASE_AF - true_af

    if calls.empty:
        return

    calls["true_af"] = true_af.apply("{:.3f}".format)
    sns.stripplot("true_af", "error", hue="caller", data=calls, palette=colors, dodge=True, jitter=True, alpha=0.5, size=2, rasterized=True)
    sns.boxplot("true_af", "error", hue="caller", data=calls, color="white", fliersize=0, linewidth=1)
    handles, labels = plt.gca().get_legend_handles_labels()
    n = len(calls.caller.unique())
    plt.legend(handles[n:], labels[n:], loc="best")

    plt.xlabel("true allele frequency")
    plt.ylabel("predicted - truth")
    plt.ylim((-1,1))
    plt.grid(axis="y")


all_calls = []
all_colors = []
for calls, caller in zip(snakemake.input.prosic_calls, snakemake.params.prosic_callers):
    label = "prosic+{}".format(caller)
    calls = pd.read_table(calls)
    calls["caller"] = label
    all_calls.append(calls)
    all_colors.append(colors[caller])

all_calls = pd.concat(all_calls)

plot(all_calls, all_colors)

plt.savefig(snakemake.output[0], bbox_inches="tight")
