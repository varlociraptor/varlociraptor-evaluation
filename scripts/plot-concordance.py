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


def common_excl(calls):
    if len(calls) == 0:
        return None
    matching = lambda c: len(c[c.MATCHING >= 0])
    run0 = calls[calls.run == snakemake.params.runs[0]]
    run1 = calls[calls.run == snakemake.params.runs[1]]

    common = matching(run0) + matching(run1)
    exclusive = len(run0) + len(run1) - common
    concordance = common / (common + exclusive)
    return common, exclusive, concordance

max_exclusive = 0

def concordance_curve(concordance):
    exclusive = np.arange(max_exclusive)
    common = concordance * exclusive / (1 - concordance)
    plt.loglog(exclusive, common, ":", color="grey")
    plt.text(max_exclusive + 10, common[-1], "{:g}".format(concordance))

concordances = []
for label, calls in all_calls.groupby("label"):
    caller = calls.caller.iloc[0]
    caller_config = snakemake.config["caller"][caller]
    is_prosic = label.startswith("prosic")
    #if is_prosic:
    #    t = -10 * np.log10(0.95)
    #    calls = calls[calls["PROB_SOMATIC"] <= t]

    #if is_prosic:
    #    thresholds = -10 * np.log10(np.array([0.75, 0.85, 0.95, 0.999])) # calls.score.quantile(np.linspace(0.0, 1.0, 50))
    #    invert = caller_config.get("invert", False)
    #    _concordances = []
    #    _commons = []
    #    for t in thresholds:
    #        if invert:
    #            c = calls[calls.score >= t]
    #        else:
    #            c = calls[calls.score <= t]
    #        common, exclusive, concordance = common_excl(c)
    #        _commons.append(common)
    #        _concordances.append(concordance)
    #
    #    plt.semilogx(_commons, _concordances, "-s", color=colors[caller], label=r"prosic+{}".format(caller))
    style = "s" if is_prosic else "o"
    common, exclusive, concordance = common_excl(calls)
    plt.semilogx(common, concordance, style, color=colors[caller], label=label)


sns.despine()
plt.xlabel("number of common calls")
plt.ylabel("concordance")
plt.legend()

plt.savefig(snakemake.output[0], bbox_inches="tight")
