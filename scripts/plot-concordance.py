from itertools import product
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import common
import numpy as np


all_calls = []

def load(call_files, runs, call_type):
    for (calls, run) in zip(call_files, runs):
        calls = common.load_variants(calls, vartype=snakemake.wildcards.vartype)
        calls["run"] = run.run
        calls["label"] = "prosic+{}".format(caller)
        calls["type"] = call_type
        all_calls.append(calls)

load(snakemake.input.prosic_calls, snakemake.params.prosic_runs, "curve")
load(snakemake.input.default_calls, snakemake.params.default_runs, "curve")
load(snakemake.input.adhoc_calls, snakemake.params.adhoc_runs, "dot")
all_calls = pd.concat(all_calls)

runs = all_calls.run.unique()

