import pandas as pd
import numpy as np


def load_variants(path, minlen, maxlen, vartype=None, constrain=None, min_af=None, max_af=None):
    variants = pd.read_table(path, header=[0, 1])

    # omit sample information
    variants = variants["VARIANT"]
    variants.index = np.arange(variants.shape[0])

    # constrain type
    if variants["SVTYPE"].isnull().any():
        if vartype == "DEL":
            variants = variants[variants["REF"].str.len() > variants["ALT"].str.len()]
        elif vartype == "INS":
            variants = variants[variants["ALT"].str.len() > variants["REF"].str.len()]
        else:
            assert False, "Unsupported variant type"
    else:
        variants = variants[variants["SVTYPE"] == vartype]

    # constrain length
    if variants["SVLEN"].isnull().any():
        if variants["END"].isnull().any():
            variants["SVLEN"] = (variants["ALT"].str.len() - variants["REF"].str.len()).abs()
            print("REF ALT comp")
        else:
            print("use END")
            variants["SVLEN"] = variants["END"] - variants["POS"]
            print(variants[["SVLEN", "POS", "END", "MATCHING"]].head())
    variants = variants[(variants["SVLEN"].abs() >= minlen) & (variants["SVLEN"].abs() < maxlen)]

    # only autosomes
    variants = variants[(variants["CHROM"] != "chrX") & (variants["CHROM"] != "chrY")]

    if constrain is not None:
        valid = (variants["MATCHING"] < 0) | (variants["MATCHING"].isin(constrain.index))
        variants = variants[valid]

    if min_af is not None and max_af is not None:
        valid = (variants["AF"] <= max_af) & (variants["AF"] >= min_af)
        variants = variants[valid]

    print("total variants", variants.shape[0])
    if "MATCHING" in variants.columns:
        print("matching variants", (variants["MATCHING"] >= 0).sum())

    return variants


def precision(calls):
    p = calls.shape[0]
    matches = calls.loc[calls["MATCHING"] >= 0, "MATCHING"]
    tp = matches.shape[0]
    precision = tp / p
    return precision


def recall(calls, truth):
    p = calls.shape[0]
    matches = calls.loc[calls.MATCHING.isin(truth.index), "MATCHING"]
    tp = matches.unique().shape[0]
    t = truth.shape[0]
    recall = tp / t
    return recall
