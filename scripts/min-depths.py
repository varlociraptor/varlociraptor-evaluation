import pandas as pd
import numpy as np


min_depth = None
for f in snakemake.input:
    d = pd.read_table(f, header=None, names=["chrom", "pos", "depth"], index_col=[0, 1], engine="c", dtype={"chrom": str, "pos": np.int32, "depth": np.int16})
    if min_depth is None:
        min_depth = d
    else:
        c = pd.concat([min_depth, d], axis=1)
        c.fillna(0, inplace=True)
        min_depth = c.min(axis=1)
        min_depth = min_depth[min_depth > 1]
min_depth.to_csv(snakemake.output[0], sep="\t")
