from itertools import product
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import common
import numpy as np


MIN_CALLS = 100

colors = common.get_colors(snakemake.config)

props = product(snakemake.params.callers,
                snakemake.params.len_ranges, snakemake.params.fdrs)

calls = []

for _calls, (caller, len_range, fdr) in zip(snakemake.input.prosic_calls, props):
    calls.append({"caller": caller, "len_range": len_range, "fdr": float(fdr), "calls": _calls})

calls = pd.DataFrame(calls)
calls = calls.set_index("caller", drop=False)


def plot_len_range(minlen, maxlen):

    def plot(caller):
        color = colors[caller]
        label = "prosic+{}".format(caller)
        fdrs = []
        alphas = []
        calls_ = calls.loc[caller]
        calls_ = calls_[calls_["len_range"].map(lambda r: r == [minlen, maxlen])]
        calls_ = calls_.sort_values("fdr")
        for e in calls_.itertuples():
            c = pd.read_table(e.calls)
            n = c.shape[0]
            if n < MIN_CALLS:
                continue
            true_fdr = 1.0 - common.precision(c)
            if fdrs and fdrs[-1] == true_fdr:
                continue
            fdrs.append(true_fdr)
            alphas.append(e.fdr)
        plt.plot(alphas, fdrs, ".-", color=color, label=label)


    for caller in calls.index:
        plot(caller)

    plt.plot([0, 1], [0, 1], ":g")

    sns.despine()
    ax = plt.gca()
    handles, _ = ax.get_legend_handles_labels()
    return ax, handles

common.plot_len_ranges(
    snakemake.params.len_ranges,
    plot_len_range,
    xlabel="FDR threshold",
    ylabel="true FDR")

plt.savefig(snakemake.output[0], bbox_inches="tight")
