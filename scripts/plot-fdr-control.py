from itertools import product
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import common
import numpy as np


MIN_CALLS = 100

colors = common.get_colors(snakemake.config)

calls = []

for _calls, (caller, fdr) in zip(snakemake.input.prosic_calls, snakemake.params.props):
    calls.append({"caller": caller, "fdr": float(fdr), "calls": _calls})

calls = pd.DataFrame(calls)
calls = calls.set_index("caller", drop=False)


def plot(caller):
    color = colors[caller]
    label = "prosic+{}".format(caller)
    fdrs = []
    alphas = []
    calls_ = calls.loc[caller].sort_values("fdr")
    for e in calls_.itertuples():
        c = pd.read_table(e.calls)
        n = c.shape[0]
        if n < MIN_CALLS:
            continue
        true_fdr = 1.0 - common.precision(c)
        if fdrs and fdrs[-1] == true_fdr:
            continue
        fdrs.append(true_fdr)
        alphas.append(e.fdr)
    plt.plot(alphas, fdrs, ".-", color=color, label=label)


for caller in calls.index:
    plot(caller)

sns.despine()
plt.legend()
plt.xlabel("FDR threshold")
plt.ylabel("true FDR")
plt.plot([0, 1], [0, 1], ":g")
plt.savefig(snakemake.output[0], bbox_inches="tight")
