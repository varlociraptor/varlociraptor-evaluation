import pandas as pd
from common import load_variants


minlen = int(snakemake.wildcards.minlen)
maxlen = int(snakemake.wildcards.maxlen)
vartype = snakemake.wildcards.vartype

truth = load_variants(snakemake.input.truth, minlen, maxlen, vartype=vartype)
calls = load_variants(snakemake.input.calls, minlen, maxlen, vartype=vartype, constrain=truth)

calls["is_tp"] = calls["MATCHING"] >= 0
calls["score"] = calls[snakemake.params.score]

calls.to_csv(snakemake.output[0], sep="\t")
