from itertools import product
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import common
import numpy as np
import math
from matplotlib.lines import Line2D

MIN_CALLS = 10

vartype = snakemake.wildcards.vartype
colors = common.get_colors(snakemake.config)


def props(callers):
    return product(callers, snakemake.params.len_ranges)


def plot_len_range(minlen, maxlen):

    truth = common.load_variants(
        snakemake.input.truth, minlen, maxlen, vartype=vartype)

    afs = pd.Series(truth.TAF.unique()).sort_values()

    def plot(calls,
             label,
             color,
             prosic=True,
             style="-.",
             markersize=4):
        calls = pd.read_table(calls, index_col=0)
        if len(calls) < 10:
            return
        if prosic:
            phred = lambda p: -10 * math.log10(p)
            def calc_recall(p):
                c = calls[calls.score <= phred(p)]
                return [common.recall(c, truth[truth.TAF >= af]) for af in afs]

            plt.fill_between(
                afs,
                calc_recall(0.98),
                calc_recall(0.8),
                color=color,
                label=label,
                alpha=0.8)
        else:
            recall = [common.recall(calls, truth[truth.TAF >= af]) for af in afs]
            plt.plot(
                afs,
                recall,
                style,
                color=color,
                label=label,
                markersize=markersize)

    handles = []
    for calls, (caller,
                len_range) in zip(snakemake.input.prosic_calls,
                                  props(snakemake.params.prosic_callers)):
        if len_range[0] != minlen and len_range[1] != maxlen:
            continue
        label = "prosic+{}".format(caller)
        plot(calls, label, colors[caller], prosic=True)
        handles.append(Line2D([0], [0], color=colors[caller], label=label))

    for calls, (caller, len_range) in zip(snakemake.input.adhoc_calls,
                             props(snakemake.params.adhoc_callers)):
        if len_range[0] != minlen and len_range[1] != maxlen:
            continue
        color = colors[caller]
        plot(calls, caller, color, style=":", prosic=False)
        handles.append(Line2D([0], [0], color=color, label=caller))

    sns.despine()
    ax = plt.gca()
    return ax, handles


common.plot_len_ranges(
    snakemake.params.len_ranges,
    plot_len_range,
    xlabel="allele frequency",
    ylabel="recall")

plt.savefig(snakemake.output[0], bbox_inches="tight")
