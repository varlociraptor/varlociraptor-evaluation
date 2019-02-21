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


class NotEnoughObservationsException(Exception):
    pass


MIN_CALLS = 10

vartype = snakemake.wildcards.vartype
colors = common.get_colors(snakemake.config)



varlociraptor_calls_low = [pd.read_table(f) for f in snakemake.input.varlociraptor_calls_low]
varlociraptor_calls_high = [pd.read_table(f) for f in snakemake.input.varlociraptor_calls_high]
adhoc_calls = [pd.read_table(f) for f in snakemake.input.adhoc_calls]



def calc_concordance(calls):
    n = len(calls)
    return (calls["concordance_count"] > 1).sum() / n


def plot_len_range(minlen, maxlen, yfunc=None, yscale=None):
    handles_varlociraptor = []
    handles_adhoc = []
    for i, caller in enumerate(snakemake.params.callers):
        def plot_calls(calls, label, color, style, calls_lower=None):
            def get_xy(calls, caseafs=None):
                svlen = calls.loc[:, calls.columns.str.startswith("SVLEN")].abs()
                # at least one of the calls has a valid svlen
                valid = ((svlen >= minlen) & (svlen <= maxlen)).sum(axis=1) >= 1
                calls = calls[valid]
                if caseafs is None:
                    caseafs = calls["max_case_af"].dropna().unique()
                y = []
                _caseafs = []
                for caseaf in sorted(caseafs):
                    _calls = calls[calls["max_case_af"] >= caseaf]
                    if len(_calls) < MIN_CALLS:
                        continue
                    _caseafs.append(caseaf)
                    y.append(yfunc(_calls))
                return _caseafs, y

            x, y = get_xy(calls)
            if not x:
                raise NotEnoughObservationsException()
            if calls_lower is not None:
                _, y2 = get_xy(calls_lower, caseafs=x)
                return plt.fill_between(x, y, y2, label=label, color=color, alpha=0.6)
            else:
                if style != "-":
                    plt.plot(x, y, "-", color="white", alpha=0.8)
                return plt.plot(x, y, style, label=label, color=color)[0]

        color = colors[snakemake.params.callers[i]]
        try:
            handles_varlociraptor.append(
                plot_calls(
                    varlociraptor_calls_high[i], 
                    "varlociraptor+{}".format(caller), 
                    color=color, style="-", 
                    calls_lower=varlociraptor_calls_low[i]))
        except NotEnoughObservationsException:
            # skip plot
            pass
        try:
            handles_adhoc.append(plot_calls(adhoc_calls[i], caller, color=color, style=":"))
        except NotEnoughObservationsException:
            # skip plot
            pass

    handles = handles_varlociraptor + handles_adhoc
    sns.despine()
    ax = plt.gca()
    if yscale is not None:
        ax.set_yscale(yscale)
    return ax, handles

plt.figure(figsize=(8, 4))
plt.subplot(121)
plot_len_range(1, 2500, yfunc=calc_concordance)
plt.xlabel("$\geq$ tumor allele frequency")
plt.ylabel("concordance")

plt.subplot(122)
ax, handles = plot_len_range(1, 2500, yfunc=lambda calls: len(calls), yscale="log")
plt.xlabel("$\geq$ tumor allele frequency")
plt.ylabel("# of calls")

ax.legend(handles=handles, loc="best")

plt.tight_layout()

plt.savefig(snakemake.output[0], bbox_inches="tight")
