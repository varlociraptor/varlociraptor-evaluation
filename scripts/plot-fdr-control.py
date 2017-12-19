from itertools import product
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import common
import numpy as np


MIN_CALLS = 100
ALPHAS = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]

minlen = int(snakemake.wildcards.minlen)
maxlen = int(snakemake.wildcards.maxlen)
vartype = snakemake.wildcards.vartype

truth = common.load_variants(snakemake.input.truth, minlen, maxlen, vartype=vartype)
colors = common.get_colors(snakemake.config)

def plot(calls, gammas, label, color, line=True, style=".-", invert=False):
    calls = pd.read_table(calls)
    gammas = pd.read_table(gammas, index_col=0, squeeze=True)
    if len(calls) < 10:
        return

    precisions = []
    alphas = []
    for alpha in ALPHAS:
        try:
            gamma = gammas[alpha]
        except KeyError:
            continue
        if gamma == 0:
            continue

        c = calls[calls.score <= gamma]
        if c.shape[0] < MIN_CALLS:
            continue

        precisions.append(common.precision(c))
        alphas.append(alpha)

    plt.plot(alphas, precisions, style, color=color, label=label)


for calls, gammas, (caller, purity) in zip(snakemake.input.prosic_calls, snakemake.input.prosic_gammas, product(snakemake.params.prosic_callers, snakemake.params.purity)):
    label = "prosic+{}".format(caller)
    if len(snakemake.params.purity) > 1:
        label += " (purity={})".format(purity)
    plot(calls, gammas, label, colors[caller])

sns.despine()
plt.legend()
plt.xlabel("FDR threshold")
plt.ylabel("precision")
plt.plot([0, 1], [1, 0], ":g")
plt.savefig(snakemake.output[0], bbox_inches="tight")
