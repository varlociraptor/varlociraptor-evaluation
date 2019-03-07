import math
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import common
import numpy as np


vartype = snakemake.wildcards.vartype
colors = common.get_colors(snakemake.config)

truth = common.load_variants(snakemake.input.truth, vartype=vartype)

all_calls = []
for caller, calls in zip(snakemake.params.callers, snakemake.input.calls):
    calls = pd.read_table(calls)
    calls.loc[:, "caller"] = caller
    all_calls.append(calls)
all_calls = pd.concat(all_calls)

def plot_depth_range(min_depth, max_depth):
    dp = all_calls["TUMOR_DP"]
    calls = all_calls[(dp >= min_depth) & (dp <= max_depth)]
    calls = calls[calls.is_tp]
    true_af = truth.loc[calls.MATCHING].reset_index().TAF
    calls = calls.reset_index()
    calls["true_af"] = true_af
    calls["error"] = calls.CASE_AF - true_af
    true_af = pd.Series(calls["true_af"].unique()).sort_values()

    # standard deviation when sampling in binomial process from allele freq
    # this is the expected sampling error within the correctly mapped fragments
    avg_depth = calls["TUMOR_DP"].mean()
    sd = true_af.apply(lambda af: 1 / avg_depth * math.sqrt(avg_depth * af * (1 - af)))
    x = np.arange(len(true_af))
    offsets = [-0.5, 0.5]
    y_upper = np.array([v for v in sd for o in offsets])
    y_lower = np.maximum(-y_upper, [-f for f in true_af for o in offsets])
    plt.fill_between([v + o for v in x for o in offsets], y_lower, y_upper, color="#EEEEEE", zorder=-5)

    calls["true_af"] = calls["true_af"].apply("{:.3f}".format)
    size = 1 if min_depth > 20 else 2
    sns.boxplot(x="true_af", y="error", data=calls, color="white", fliersize=0, linewidth=1)
    sns.stripplot("true_af", "error", hue="caller", data=calls, palette=colors, dodge=True, jitter=True, alpha=0.5, size=size, rasterized=True)
    #sns.scatterplot(x="true_af", y="CASE_AF", hue="caller", hue_order=colors.keys(), palette=colors.values(), data=calls, rasterized=True)
    
    #sd = lambda af, n: 1 / n * math.sqrt(n * af * (1 - af))
    #all_sd = lambda n: np.array([sd(af, n) for af in true_af])
    # lower
    #plt.fill_between(true_af, true_af - all_sd(avg_depth), true_af + all_sd(avg_depth), color="grey")

    sns.despine()
    plt.xticks(rotation="vertical")
    ax = plt.gca()
    ax.legend().remove()
    handles, labels = ax.get_legend_handles_labels()

    return ax, handles

common.plot_ranges(
    snakemake.params.depth_ranges,
    plot_depth_range,
    "true allele frequency",
    "predicted - truth")

plt.savefig(snakemake.output[0], bbox_inches="tight")
