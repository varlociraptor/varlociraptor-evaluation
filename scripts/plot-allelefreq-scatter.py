import math
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import common
import numpy as np

MIN_COUNT=20

vartype = snakemake.wildcards.vartype
colors = common.get_colors(snakemake.config)

truth = common.load_variants(snakemake.input.truth, vartype=vartype)

all_calls = []
for caller, calls in zip(snakemake.params.callers, snakemake.input.calls):
    calls = pd.read_table(calls)
    calls.loc[:, "caller"] = caller
    all_calls.append(calls)
all_calls = pd.concat(all_calls)

def plot(af, _):
    constrain_lower = lambda error: np.maximum(error, -af)
    constrain_upper = lambda error: np.minimum(error, 1.0 - af)

    dp = all_calls["TUMOR_DP"]
    calls = all_calls[all_calls.is_tp]
    true_af = truth.loc[calls.MATCHING].reset_index().TAF
    calls = calls.reset_index()
    calls["true_af"] = true_af
    calls = calls[calls["true_af"] == af]
    calls["error"] = calls.CASE_AF - true_af

    sns.kdeplot(calls["TUMOR_DP"], calls["error"], cmap="Blues", n_levels=50, shade=True, alpha=0.7, shade_lowest=False) #alpha=0.5, clip=((0.0, 1.0), (0.0, af)))
    plt.plot(calls["TUMOR_DP"], calls["error"], ",", color="k", lw=0, alpha=1.0, rasterized=True)
    by_depth = calls.groupby("TUMOR_DP")["error"].describe().reset_index()
    by_depth["-std"] = constrain_lower(-by_depth["std"])
    by_depth["std"] = constrain_upper(by_depth["std"])
    by_depth = by_depth[by_depth["count"] >= MIN_COUNT]
    plt.plot(by_depth.TUMOR_DP, by_depth["std"], "--", color="k")
    plt.plot(by_depth.TUMOR_DP, by_depth["-std"], "--", color="k")
    plt.plot(by_depth.TUMOR_DP, by_depth["mean"], "-", color="k")

    max_depth = calls["TUMOR_DP"].max()
    depths = np.arange(0, max_depth)
    # standard deviation when sampling in binomial process from allele freq
    # this is the expected sampling error within the correctly mapped fragments
    sd = np.array([1.0 / depth * math.sqrt(depth * af * (1.0 - af)) for depth in depths])
    plt.fill_between(depths, constrain_lower(-sd), constrain_upper(sd), color="grey", alpha=0.5)
    
    sns.despine()
    plt.xticks(rotation="vertical")
    ax = plt.gca()
    ax.legend().remove()
    handles, labels = ax.get_legend_handles_labels()
    plt.ylim((-1.0, 1.0))
    plt.xlim((0, 60))

    return ax, []

afs = [(af, af) for af in truth.TAF.sort_values().unique()]

common.plot_ranges(
    afs,
    plot,
    "depth",
    "predicted - truth")

plt.savefig(snakemake.output[0], bbox_inches="tight")
