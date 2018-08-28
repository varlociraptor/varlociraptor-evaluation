import pandas as pd
import numpy as np
from common import load_variants


minlen = int(snakemake.wildcards.minlen)
maxlen = int(snakemake.wildcards.maxlen)
vartype = snakemake.wildcards.vartype

if snakemake.wildcards.mode == "prosic":
    score = snakemake.config["caller"]["prosic"]["score"]
    # calls are already filtered by FDR control step
    minlen = None
    maxlen = None
elif snakemake.wildcards.mode == "default":
    score = snakemake.config["caller"][snakemake.wildcards.caller]["score"]
else:
    score = None

calls = load_variants(snakemake.input.calls, vartype=vartype, minlen=minlen, maxlen=maxlen)

calls["is_tp"] = calls["MATCHING"] >= 0

calls["score"] = calls[score] if score else np.nan

calls.to_csv(snakemake.output[0], sep="\t")
