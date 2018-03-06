from itertools import product
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import common
import numpy as np


MIN_CALLS = 10


minlen = int(snakemake.wildcards.minlen)
maxlen = int(snakemake.wildcards.maxlen)
vartype = snakemake.wildcards.vartype

truth = common.load_variants(snakemake.input.truth, minlen, maxlen, vartype=vartype)
colors = common.get_colors(snakemake.config)


def plot(calls, label, color, line=True, style="-", invert=False):
    calls = pd.read_table(calls)
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
    else:
        precision = [common.precision(calls)]
        recall = [common.recall(calls, truth)]
        style = "."
        print(label, calls.shape[0], precision, recall)
    plt.plot(recall, precision, style, color=color, label=label)


for calls, (caller, purity) in zip(snakemake.input.prosic_calls, product(snakemake.params.prosic_callers, snakemake.params.purity)):
    label = "prosic+{}".format(caller)
    if len(snakemake.params.purity) > 1:
        label += " (purity={})".format(purity)
    plot(calls, label, colors[caller])

for calls, caller in zip(snakemake.input.default_calls, snakemake.params.default_callers):
    plot(calls, caller, colors[caller], style=":", invert=snakemake.config["caller"][caller].get("invert", False))

for calls, caller in zip(snakemake.input.adhoc_calls, snakemake.params.adhoc_callers):
    plot(calls, caller, colors[caller], line=False)

sns.despine()
plt.legend(loc="lower right")
plt.xlabel("recall")
plt.ylabel("precision")
plt.savefig(snakemake.output[0], bbox_inches="tight")
