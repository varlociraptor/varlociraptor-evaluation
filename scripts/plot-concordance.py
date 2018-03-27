from itertools import product
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import common
import numpy as np


all_calls = []

def load(call_files, runs, call_type, label):
    for (calls, run) in zip(call_files, runs):
        calls = common.load_variants(calls, vartype=snakemake.wildcards.vartype)
        calls["run"] = run.run
        calls["label"] = label(caller)
        calls["type"] = call_type
        all_calls.append(calls)

load(snakemake.input.prosic_calls, snakemake.params.prosic_runs, "curve", "prosic+{}".format)
load(snakemake.input.default_calls, snakemake.params.default_runs, "curve", lambda caller: caller)
load(snakemake.input.adhoc_calls, snakemake.params.adhoc_runs, "dot", "{}-adhoc".format)
all_calls = pd.concat(all_calls)

all_calls.groupby("label") # TODO go on
