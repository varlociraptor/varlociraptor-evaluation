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

    def plot(calls,
             label,
             color,
             line=True,
             style="-",
             invert=False,
             markersize=4):
        calls = pd.read_table(calls, index_col=0)
        if len(calls) < 10:
            return
        if line:
            thresholds = calls.score.quantile(np.linspace(0.0, 1.0, 50))
            precision = []
            recall = []
            for t in thresholds:
                if invert:
                    c = calls[calls.score >= t]
                else:
                    c = calls[calls.score <= t]
                p = common.precision(c)
                r = common.recall(c, truth)
                print(label, t, c.shape[0], p, r)
                if len(c) < 10:
                    print(c)
                precision.append(p)
                recall.append(r)
            if len(precision) <= 2:
                print("skipping curve because we have too few values")
                return
        else:
            precision = [common.precision(calls)]
            recall = [common.recall(calls, truth)]
            style = "."
            print(label, calls.shape[0], precision, recall)
        plt.plot(
            recall,
            precision,
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
        plot(calls, label, colors[caller])
        handles.append(Line2D([0], [0], color=colors[caller], label=label))

    for calls, (caller,
                len_range) in zip(snakemake.input.default_calls,
                                  props(snakemake.params.default_callers)):
        if len_range[0] != minlen and len_range[1] != maxlen:
            continue
        color = colors[caller]
        plot(
            calls,
            caller,
            color,
            style=":",
            invert=snakemake.config["caller"][caller].get("invert", False))
        if caller in snakemake.params.adhoc_callers:
            handles.append(Line2D([0], [0], markersize=10, markerfacecolor=color, markeredgecolor=color, color=color, label=caller, marker=".", linestyle=":"))
        else:
            handles.append(Line2D([0], [0], color=color, label=caller, linestyle=":"))

    for calls, (caller, len_range) in zip(snakemake.input.adhoc_calls,
                             props(snakemake.params.adhoc_callers)):
        if len_range[0] != minlen and len_range[1] != maxlen:
            continue
        color = colors[caller]
        plot(calls, caller, color, markersize=10, line=False)
        if caller not in snakemake.params.default_callers:
            handles.append(Line2D([0], [0], markersize=10, markerfacecolor=color, markeredgecolor=color, label=caller, marker=".", lw=0))

    sns.despine()
    ax = plt.gca()
    return ax, handles


common.plot_len_ranges(
    snakemake.params.len_ranges,
    plot_len_range,
    xlabel="recall",
    ylabel="precision")

plt.savefig(snakemake.output[0], bbox_inches="tight")
