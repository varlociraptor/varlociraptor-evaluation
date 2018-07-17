from itertools import product
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import common
import numpy as np

minlen = int(snakemake.wildcards.minlen)
maxlen = int(snakemake.wildcards.maxlen)
colors = common.get_colors(snakemake.config)


all_calls = []

def load(call_files, runs, call_type, label):
    for (calls, (caller, run)) in zip(call_files, runs):
        calls = common.load_variants(calls, vartype=snakemake.wildcards.vartype, minlen=minlen, maxlen=maxlen)
        if calls.empty:
            continue

        calls["run"] = run
        calls["label"] = label(caller)
        calls["caller"] = caller
        caller_config = snakemake.config["caller"][caller]
        if "score" in caller_config:
            calls["score"] = calls[caller_config["score"]]
        calls["type"] = call_type
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
    if len(run0) < 10 and len(run1) < 10:
        return None
    common = min(matching(run0), matching(run1))
    print(calls.caller.iloc[0], calls.shape, len(run0), len(run1), common)
    return common / (len(run0) + len(run1) - common)


for label, calls in all_calls.groupby("label"):
    caller = calls.caller.iloc[0]
    caller_config = snakemake.config["caller"][caller]
    is_prosic = label.startswith("prosic")

    if is_prosic or "score" in caller_config:
        thresholds = calls.score.quantile(np.linspace(0.0, 1.0, 50))
        invert = caller_config.get("invert", False)
        concordance = []
        counts = []
        for t in thresholds:
            if invert:
                c = calls[calls.score >= t]
            else:
                c = calls[calls.score <= t]
            j = jaccard(c)
            if j is not None:
                concordance.append(j)
                counts.append(len(c))
        s = "-" if is_prosic else ":"
        plt.semilogx(counts, concordance, s, label=label, color=colors[caller])
    if not is_prosic and caller_config.get("adhoc", False):
        j = jaccard(calls)
        if j is not None:
            plt.semilogx(len(calls), j, "o", color=colors[caller])

sns.despine()
plt.xlabel("number of calls")
plt.ylabel("concordance")
plt.legend()

plt.savefig(snakemake.output[0], bbox_inches="tight")
