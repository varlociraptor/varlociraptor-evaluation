from itertools import product
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import common
import numpy as np
import math
from matplotlib.lines import Line2D

MIN_CALLS = 10

vartype = snakemake.wildcards.vartype
colors = common.get_colors(snakemake.config)



prosic_calls = [pd.read_table(f) for f in snakemake.input.prosic_calls]
adhoc_calls = [pd.read_table(f) for f in snakemake.input.adhoc_calls]



def calc_concordance(calls):
    n = len(calls)
    return (calls["concordance_count"] > 1).sum() / n


def plot_len_range(minlen, maxlen, yfunc=None, pltfunc=plt.plot):
    for i, caller in enumerate(snakemake.params.callers):
        def plot_calls(calls, label, color, style):
            svlen = calls.loc[:, calls.columns.str.startswith("SVLEN")].abs()
            # at least one of the calls has a valid svlen
            valid = ((svlen >= minlen) & (svlen <= maxlen)).sum(axis=1) >= 1
            calls = calls[valid]
            caseafs = calls["max_case_af"].unique()
            y = []
            n_calls = []
            _caseafs = []
            for caseaf in sorted(caseafs):
                _calls = calls[calls["max_case_af"] >= caseaf]
                if len(_calls) < MIN_CALLS:
                    continue
                _caseafs.append(caseaf)
                y.append(yfunc(_calls))
                n = len(_calls)
                n_calls.append(n)

            pltfunc(_caseafs, y, style, label=label, color=color)

        color = colors[snakemake.params.callers[i]]
        plot_calls(prosic_calls[i], "prosic+{}".format(caller), color=color, style="-")
        plot_calls(adhoc_calls[i], caller, color=color, style=":")
            
    sns.despine()
    ax = plt.gca()
    handles, _ = ax.get_legend_handles_labels()
    return ax, handles

plt.figure(figsize=(4, 7))
plt.subplot(211)
plot_len_range(1, 2500, yfunc=calc_concordance)
plt.ylabel("concordance")

plt.subplot(212)
plot_len_range(1, 2500, yfunc=lambda calls: len(calls), pltfunc=plt.semilogy)
plt.xlabel("tumor allele frequency")
plt.ylabel("# of calls")

plt.tight_layout()

plt.savefig(snakemake.output[0], bbox_inches="tight")
