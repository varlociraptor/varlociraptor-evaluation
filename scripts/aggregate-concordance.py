from common import load_variants
import networkx as nx
import pandas as pd
import numpy as np

vartype = snakemake.wildcards.vartype

all_variants = [load_variants(f, vartype=vartype) for f in snakemake.input.calls]

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

if snakemake.wildcards.mode != "varlociraptor":
    varlociraptor_variants = [load_variants(f, vartype=vartype) for f in snakemake.input.varlociraptor_calls]
    for calls in varlociraptor_variants:
        calls.set_index(["CHROM", "POS", "REF", "ALT", "SVLEN"], inplace=True)
    varlociraptor_representatives = {snakemake.params.dataset_combinations[i][0]: calls for i, calls in enumerate(varlociraptor_variants)}

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
is_varlociraptor = False
for dataset_id, calls in representatives.items():
    cols = ["CHROM", "POS", "REF", "ALT", "SVLEN"]
    if "CASE_AF" in calls.columns:
        cols.extend(["CASE_AF", "PROB_SOMATIC_TUMOR"])
        is_varlociraptor = True
    calls = calls[cols]
    if snakemake.wildcards.mode != "varlociraptor":
        caseaf = calls.set_index(cols, drop=False).join(varlociraptor_representatives[dataset_id][["CASE_AF"]], how="left")["CASE_AF"]
        caseaf = caseaf[~caseaf.index.duplicated()]
        calls["CASE_AF"] = caseaf.values

    calls.columns = [c + suffix(dataset_name(dataset_id)) for c in calls.columns]
    if aggregated is None:
        aggregated = calls
    else:
        aggregated = aggregated.join(calls, how="outer", lsuffix="", rsuffix="")

# Forget the component id. Otherwise, we might run into errors with duplicate elements 
# in the index below. These can occur if there are multiple ambiguous calls.
aggregated.reset_index(inplace=True, drop=True)

pos_cols = aggregated.columns[aggregated.columns.str.startswith("POS_")]
is_called = (~aggregated[pos_cols].isnull()).astype(int)
is_called.columns = pos_cols.str.replace("POS_", "")
aggregated = aggregated.join(is_called, lsuffix="", rsuffix="")

aggregated.insert(len(aggregated.columns), "concordance_count", is_called.sum(axis=1))

aggregated["max_case_af"] = aggregated[aggregated.columns[aggregated.columns.str.startswith("CASE_AF")]].max(axis=1)
if is_varlociraptor:
    aggregated["max_prob_somatic_tumor"] =  aggregated[aggregated.columns[aggregated.columns.str.startswith("PROB_SOMATIC")]].min(axis=1)

aggregated.to_csv(snakemake.output[0], sep="\t", index=False)
