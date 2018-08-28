from itertools import product
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import common
import numpy as np


MIN_CALLS = 10

vartype = snakemake.wildcards.vartype
colors = common.get_colors(snakemake.config)

def props(callers):
    return product(callers, snakemake.params.len_ranges)

def plot_len_range(minlen, maxlen):

    truth = common.load_variants(snakemake.input.truth, minlen, maxlen, vartype=vartype)

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

        plt.ylim((-1,1))
        plt.grid(axis="y", linestyle=":", color="g")
        plt.despine()

        return plt.gca(), handles[n:]

    all_calls = []
    all_colors = []
    for calls, (caller, len_range) in zip(snakemake.input.prosic_calls, props(snakemake.params.prosic_callers)):
        if len_range[0] != minlen and len_range[1] != maxlen:
            continue
        label = "prosic+{}".format(caller)
        calls = pd.read_table(calls)
        calls["caller"] = label
        if not calls.empty:
            all_calls.append(calls)
            all_colors.append(colors[caller])

    all_calls = pd.concat(all_calls)

    plot(all_calls, all_colors)

common.plot_len_ranges(
    snakemake.params.len_ranges,
    plot_len_range,
    xlabel="true allele frequency",
    ylabel="predicted - truth")

plt.savefig(snakemake.output[0], bbox_inches="tight")
