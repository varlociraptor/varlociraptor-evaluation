from itertools import product
import math
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

truth = common.load_variants(snakemake.input.truth, vartype=vartype)

def props(callers):
    return product(callers, snakemake.params.len_ranges)

def plot_len_range(minlen, maxlen):


    def plot(calls, colors):
        calls = calls[calls.is_tp]
        true_af = truth.loc[calls.MATCHING].reset_index().TAF
        calls = calls.reset_index()
        calls["error"] = calls.CASE_AF - true_af

        if calls.empty:
            return

        calls["true_af"] = true_af
        true_af = pd.Series(calls["true_af"].unique()).sort_values()
        # standard deviation when sampling in binomial process from allele freq
        # this is the expected sampling error within the correctly mapped fragments
        sd = true_af.apply(lambda af: math.sqrt(af * (1 - af)))
        x = np.arange(len(true_af))
        offsets = [-0.5, 0.5]
        y_upper = np.array([v for v in sd for o in offsets])
        y_lower = np.maximum(-y_upper, [-f for f in true_af for o in offsets])
        plt.fill_between([v + o for v in x for o in offsets], y_lower, y_upper, color="#EEEEEE", zorder=-5)

        calls["true_af"] = calls["true_af"].apply("{:.3f}".format)

        size = 1 if maxlen == 30 else 2
        sns.stripplot("true_af", "error", hue="caller", data=calls, palette=colors, dodge=True, jitter=True, alpha=0.5, size=size, rasterized=True)
        sns.boxplot("true_af", "error", hue="caller", data=calls, color="white", fliersize=0, linewidth=1)

        handles, labels = plt.gca().get_legend_handles_labels()
        n = len(calls.caller.unique())

        plt.ylim((-1,1))
        plt.grid(axis="y", linestyle=":", color="grey")
        sns.despine()
        plt.xticks(rotation="vertical")
        ax = plt.gca()
        ax.legend().remove()

        return ax, handles[n:]

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

    return plot(all_calls, all_colors)

common.plot_len_ranges(
    snakemake.params.len_ranges,
    plot_len_range,
    xlabel="true allele frequency",
    ylabel="predicted - truth")

plt.savefig(snakemake.output[0], bbox_inches="tight")
