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
vartype = snakemake.wildcards.vartype

truth = common.load_variants(snakemake.input.truth, minlen, maxlen, vartype=vartype)
colors = common.get_colors(snakemake.config)

def plot(calls, label, color, line=True, style="-"):
    calls = pd.read_table(calls)
    if line:
        thresholds = calls.score.quantile(np.linspace(0.0, 1.0, 50))
        precision = []
        recall = []
        for t in thresholds:
            c = calls[calls.score <= t]
            p = common.precision(c)
            r = common.recall(c, truth)
            precision.append(p)
            recall.append(r)
    else:
        precision = [common.precision(calls)]
        recall = [common.recall(calls, truth)]
        style = "."
    plt.plot(recall, precision, style, color=color, label=label)


for calls, (caller, purity) in zip(snakemake.input.prosic_calls, product(snakemake.params.prosic_callers, snakemake.params.purity)):
    label = "prosic+{}".format(caller)
    if len(snakemake.params.purity) > 1:
        label += " (purity={})".format(purity)
    plot(calls, label, colors[caller])

for calls, caller in zip(snakemake.input.default_calls, snakemake.params.default_callers):
    plot(calls, caller, colors[caller], style=":")

for calls, caller in zip(snakemake.input.adhoc_calls, snakemake.params.adhoc_callers):
    plot(calls, caller, colors[caller], line=False)

sns.despine()
plt.legend()
plt.xlabel("recall")
plt.ylabel("precision")
plt.savefig(snakemake.output[0], bbox_inches="tight")
