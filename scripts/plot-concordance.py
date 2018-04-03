from itertools import product
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import common
import numpy as np

print(snakemake.params.len_ranges)
len_ranges = pd.Series(["{} - {}".format(*i) for i in snakemake.params.len_ranges],
                       index = pd.IntervalIndex.from_intervals([pd.Interval(*i, closed="left") for i in snakemake.params.len_ranges]))
colors = common.get_colors(snakemake.config)

all_calls = []

def load(call_files, runs, call_type, label):
    for (calls, (caller, run)) in zip(call_files, runs):
        calls = common.load_variants(calls, vartype=snakemake.wildcards.vartype)
        if calls.empty:
            continue
        calls["run"] = run
        calls["label"] = label(caller)
        calls["caller"] = caller
        calls["type"] = call_type
        len_range = len_ranges[calls["SVLEN"]]
        #calls = calls[~pd.isnull(len_range.index)]
        #calls["len"] = len_ranges[calls["SVLEN"]].values
        calls["len"] = calls["SVLEN"]
        all_calls.append(calls)

load(snakemake.input.prosic_calls, snakemake.params.prosic_runs, "prosic", "prosic+{}".format)
load(snakemake.input.adhoc_calls, snakemake.params.adhoc_runs, "adhoc", "{}".format)
all_calls = pd.concat(all_calls)


def jaccard(calls):
    if len(calls) == 0:
        return None
    matching = lambda c: len(c[c.MATCHING >= 0])
    run0 = calls[calls.run == snakemake.params.runs[0]]
    run1 = calls[calls.run == snakemake.params.runs[1]]
    if len(run0) < 50 and len(run1) < 50:
        return None
    common = min(matching(run0), matching(run1))
    return common / (len(run0) + len(run1) - common)


for label, calls in all_calls.groupby("label"):
    caller = calls.caller.iloc[0]
    def len_vs_jaccard(calls):
        lens = calls.len.sort_values().unique()
        d = pd.DataFrame({"jaccard": [jaccard(calls[calls.len <= l]) for l in lens],
                          "len": lens})
        d = d[~pd.isnull(d.jaccard)]
        return d
    if label.startswith("prosic"):
        thresholds = [common.phred_scale(p) for p in [0.999999999, 0.9999, 0.999, 0.99, 0.95]]
        #thresholds = calls.PROB_SOMATIC.quantile(np.linspace(0.0, 1.0, 5))
        for t in thresholds:
            c = calls[calls.PROB_SOMATIC <= t]
            d = len_vs_jaccard(c)
            plt.plot("len", "jaccard", "-", data=d, label=label, color=colors[caller])
    else:
        d = len_vs_jaccard(calls)
        plt.plot("len", "jaccard", "--", data=d, label=label, color=colors[caller])
sns.despine()
plt.xlabel("length")
plt.ylabel("concordance")
plt.legend()

plt.savefig(snakemake.output[0], bbox_inches="tight")
