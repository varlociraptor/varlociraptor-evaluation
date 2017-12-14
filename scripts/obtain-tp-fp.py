import pandas as pd

calls = pd.read_table(snakemake.input.calls)
truth = pd.read_table(snakemake.input.truth)
minlen = int(snakemake.wildcards.minlen)
maxlen = int(snakemake.wildcards.maxlen)
vartype = snakemake.wildcards.vartype


