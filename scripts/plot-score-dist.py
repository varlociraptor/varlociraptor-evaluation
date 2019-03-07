from itertools import product
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import common
import numpy as np
import math

vartype = snakemake.wildcards.vartype
colors = common.get_colors(snakemake.config)

def props(callers):
    return product(callers, snakemake.params.len_ranges)

phred_to_log_factor = -0.23025850929940456
log_to_phred_factor = -4.3429448190325175

def plot_len_range(minlen, maxlen):
    for calls, (caller, len_range) in zip(snakemake.input.varlociraptor_calls, props(snakemake.params.varlociraptor_callers)):
        if len_range[0] != minlen and len_range[1] != maxlen:
            continue
        label = "varlociraptor+{}".format(caller)
        calls = pd.read_table(calls)
        calls["caller"] = label
        if not calls.empty:
            color = colors[caller]
            sns.kdeplot(calls[calls.is_tp].PROB_SOMATIC_TUMOR.map(np.log), color=color, label=label)
            sns.kdeplot(calls[~calls.is_tp].PROB_SOMATIC_TUMOR.map(np.log), color=color, linestyle=":", label="")

    ax = plt.gca()
    fmt_ticks = lambda ticks: ["{:.1g}".format(np.exp(t)) for t in ticks]
    ax.set_xticklabels(fmt_ticks(plt.xticks()[0]))
    ax.legend().remove()
    handles, _ = ax.get_legend_handles_labels()
    sns.despine()

    return ax, handles

common.plot_ranges(
    snakemake.params.len_ranges,
    plot_len_range,
    xlabel="Pr(somatic) (PHRED)",
    ylabel="density")

plt.savefig(snakemake.output[0], bbox_inches="tight")
