import pandas as pd
import numpy as np
from common import load_variants


minlen = int(snakemake.wildcards.minlen)
maxlen = int(snakemake.wildcards.maxlen)
vartype = snakemake.wildcards.vartype

truth = load_variants(snakemake.input.truth, minlen, maxlen, vartype=vartype)
calls = load_variants(snakemake.input.calls, minlen, maxlen, vartype=vartype, constrain=truth)

calls["is_tp"] = calls["MATCHING"] >= 0

if snakemake.wildcards.mode == "prosic":
    score = snakemake.config["caller"]["prosic"]["score"]
elif snakemake.wildcards.mode == "default":
    score = snakemake.config["caller"][snakemake.wildcards.caller]["score"]
else:
    score = None

calls["score"] = calls[score] if score else np.nan

calls.to_csv(snakemake.output[0], sep="\t")
