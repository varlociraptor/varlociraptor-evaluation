from itertools import product
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import common


vartype = snakemake.wildcards.vartype
colors = common.get_colors(snakemake.config)

def props(callers):
    return product(callers, snakemake.params.len_ranges)

def plot_len_range(minlen, maxlen):
    for calls, (caller, len_range) in zip(snakemake.input.prosic_calls, props(snakemake.params.prosic_callers)):
        if len_range[0] != minlen and len_range[1] != maxlen:
            continue
        label = "prosic+{}".format(caller)
        calls = pd.read_table(calls)
        if not calls.empty:
            color = colors[caller]
            sns.kdeplot(calls[calls.is_tp].PROB_SOMATIC, color=color, label=label)
            sns.kdeplot(calls[~calls.is_tp]).PROB_SOMATIC, color=color, linestyle="--")
    ax = plt.gca()
    handles, _ = ax.get_legend_handles_labels()

    return ax, handles

common.plot_len_ranges(
    snakemake.params.len_ranges,
    plot_len_range,
    xlabel="posterior probability",
    ylabel="density")

plt.savefig(snakemake.output[0], bbox_inches="tight")
