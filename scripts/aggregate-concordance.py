from common import load_variants
import networkx as nx
import pandas as pd
import numpy as np

vartype = snakemake.wildcards.vartype

all_variants = [load_variants(f, vartype=vartype) for f in snakemake.input]

G = nx.Graph()
for calls, (i, j) in zip(all_variants, snakemake.params.dataset_combinations):
    calls["component"] = None
    for call in calls.itertuples():
        a = (i, call.Index)
        G.add_node(a)
        if call.MATCHING >= 0:
            b = (j, call.MATCHING)
            G.add_node(b)
            G.add_edge(a, b)

# get a set of calls for each dataset (we don't need all pairwise comparisons for that)
representatives = {snakemake.params.dataset_combinations[i][0]: calls for i, calls in enumerate(all_variants)}

# annotate calls with their component, i.e. their equivalence class
for component_id, component in enumerate(nx.connected_components(G)):
    for i, k in component:
        representatives[i].loc[k, "component"] = component_id
for calls in representatives.values():
    calls["component"] = calls["component"].astype(np.float32)
    calls.set_index("component", inplace=True)

# join calls based on their equivalence class
aggregated = None
suffix = "_{}".format
dataset_name = lambda i: snakemake.params.datasets[i]
is_prosic = False
for dataset_id, calls in representatives.items():
    cols = ["CHROM", "POS", "REF", "ALT", "SVLEN"]
    if "CASE_AF" in calls.columns:
        cols.extend(["CASE_AF", "PROB_SOMATIC_TUMOR"])
        is_prosic = True
    calls = calls[cols]
    calls.columns = [c + suffix(dataset_name(dataset_id)) for c in calls.columns]
    if aggregated is None:
        aggregated = calls
    else:
        aggregated = aggregated.join(calls, how="outer", lsuffix="", rsuffix="")

pos_cols = aggregated.columns[aggregated.columns.str.startswith("POS_")]
is_called = (~aggregated[pos_cols].isnull()).astype(int)
is_called.columns = pos_cols.str.replace("POS_", "")
aggregated = aggregated.join(is_called, lsuffix="", rsuffix="")

aggregated["concordance_count"] = is_called.sum(axis=1)
if is_prosic:
    aggregated["max_case_af"] = aggregated[aggregated.columns[aggregated.columns.str.startswith("CASE_AF")]].max(axis=1)
    aggregated["max_prob_somatic_tumor"] =  aggregated[aggregated.columns[aggregated.columns.str.startswith("PROB_SOMATIC")]].min(axis=1)

aggregated.to_csv(snakemake.output[0], sep="\t", index=False)
