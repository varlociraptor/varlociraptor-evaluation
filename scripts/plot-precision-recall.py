from itertools import product

import pandas as pd
import common


minlen = int(snakemake.wildcards.minlen)
maxlen = int(snakemake.wildcards.maxlen)
vartype = snakemake.wildcards.vartype

truth = common.load_variants(snakemake.input.truth, minlen, maxlen, vartype=vartype)
colors = common.get_colors(config)

def plot(calls, label, color, line=True, style="-"):
    if line:
        thresholds = np.linspace(calls.score.min(), calls.score.max(), 50)
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
    plt.plot(recall, precision, style, color=color)


for calls, (caller, purity) in zip(snakemake.input.prosic_calls, product(snakemake.params.prosic_callers, snakemake.params.purity)):
    label = "prosic+{}".format(caller)
    if len(snakemake.params.purity) > 1:
        label += " (purity={})".format(purity)
    plot(calls, label, common.colors[caller])

for calls, caller in zip(snakemake.input.default_calls, snakemake.params.default_callers):
    plot(calls, caller, common.colors[caller], style=":")

for calls, caller in zip(snakemake.input.adhoc_calls, snakemake.params.default_callers):
    plot(calls, caller, common.colors[caller], line=False)

sns.despine()
plt.legend()
plt.xlabel("recall")
plt.ylabel("precision")
plt.savefig(snakemake.output[0], bbox_inches="tight")
