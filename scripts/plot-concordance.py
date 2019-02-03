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



def plot_len_range(minlen, maxlen):
    for i, caller in enumerate(snakemake.params.callers):
        def plot_calls(calls, label, color, style):
            svlen = calls.loc[:, calls.columns.str.startswith("SVLEN")].abs()
            # at least one of the calls has a valid svlen
            valid = ((svlen >= minlen) & (svlen <= maxlen)).sum(axis=1) >= 1
            calls = calls[valid]
            caseafs = calls["max_case_af"].unique()
            concordance = []
            n_calls = []
            _caseafs = []
            for caseaf in sorted(caseafs):
                _calls = calls[calls["max_case_af"] >= caseaf]
                if len(_calls) < MIN_CALLS:
                    continue
                _caseafs.append(caseaf)
                n = len(_calls)
                concordance.append((_calls["concordance_count"] > 1).sum() / n)
                n_calls.append(n)

            plt.plot(_caseafs, concordance, style, label=label, color=color)

        color = colors[snakemake.params.callers[i]]
        plot_calls(prosic_calls[i], "prosic+{}".format(caller), color=color, style="-")
        plot_calls(adhoc_calls[i], caller, color=color, style=":")
            
    sns.despine()
    ax = plt.gca()
    handles, _ = ax.get_legend_handles_labels()
    return ax, handles


common.plot_len_ranges(
    snakemake.params.len_ranges,
    plot_len_range,
    xlabel="tumor allele frequency",
    ylabel="concordance")

plt.savefig(snakemake.output[0], bbox_inches="tight")
