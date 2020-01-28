import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import numpy as np
import seaborn as sns
import math


def plot_ranges(ranges, plot_range, xlabel, ylabel, row_offset=0, nrow_offset=0, row_span=1, fig=None, gs=None, legend=True, fig_height=None, legend_outside=False):
    ncols = 4 if len(ranges) == 4 else min(3, len(ranges))
    nrows = int(math.ceil(len(ranges) / ncols)) + nrow_offset
    if fig is None:
        if fig_height is None:
            fig_height = 4 * nrows
        fig_width = 4 * ncols
        if legend_outside:
            fig_width += 2.5
        fig = plt.figure(figsize=(fig_width, fig_height))
        gs = gridspec.GridSpec(nrows, ncols, figure=fig)
    axes = []
    all_handles = []
    seen = set()
    for i, (lower, upper) in enumerate(ranges):
        row = i // ncols + row_offset
        col = i % ncols
        fig.add_subplot(gs[row:row+row_span, col]) 
        ax, handles = plot_range(lower, upper)

        if col == 0:
            plt.ylabel(ylabel)
        else:
            plt.ylabel("")
        if row + row_span == nrows:
            plt.xlabel(xlabel)
        else:
            plt.xlabel("")

        if row_offset == 0 and len(ranges) > 1:
            if lower == upper:
                if isinstance(lower, float):
                    lower = "{:.3g}".format(lower)
                plt.title(lower)
            else:
                plt.title("{} - {}".format(lower, upper))

        axes.append(ax)
        for handle in handles:
            label = handle.get_label()
            if label not in seen:
                seen.add(label)
                all_handles.append(handle)

    if legend:
        if legend_outside:
            axes[-1].legend(handles=all_handles, loc="upper left", bbox_to_anchor=(1.01, 1.0))
        else:
            axes[0].legend(handles=all_handles, loc="best")
    plt.tight_layout()
    return fig, gs


def load_variants(path,
                  minlen=None,
                  maxlen=None,
                  vartype=None,
                  constrain=None,
                  min_af=None,
                  max_af=None):
    variants = pd.read_table(path, header=[0, 1])

    # store tumor AF estimate in CASE_AF column
    try:
        case_af = variants.loc[:, ("tumor", "AF")]
        variants.loc[:, ("VARIANT", "CASE_AF")] = case_af
    except KeyError:
        # ignore if no AF estimate for tumor is present
        pass
    try:
        dp = variants.loc[:, ("tumor", "DP")]
        variants.loc[:, ("VARIANT", "TUMOR_DP")] = dp
    except KeyError:
        # ignore if not present
        pass

    variants = variants["VARIANT"]
    variants["CHROM"] = variants["CHROM"].astype(str)

    variants.index = np.arange(variants.shape[0])

    # constrain type
    if vartype == "DEL":
        is_allele_del = (variants["REF"].str.len() > 1) & (variants["ALT"].str.len() == 1)
        is_sv_del = variants["ALT"] == "<DEL>"
        isdel = is_allele_del | is_sv_del

        if "SVTYPE" in variants.columns:
            variants = variants[(variants["SVTYPE"].astype(str) == "DEL")
                                | (isdel & variants["SVTYPE"].isnull())]
        else:
            variants = variants[isdel]
    elif vartype == "INS":
        isins = (variants["REF"].str.len() == 1) & (variants["ALT"].str.len() >
                                                    1)
        if "SVTYPE" in variants.columns:
            variants = variants[(variants["SVTYPE"].astype(str) == "INS")
                                | (isins & variants["SVTYPE"].isnull())]
        else:
            variants = variants[isins]
    else:
        assert False, "Unsupported variant type"

    # constrain length
    if "SVLEN" not in variants.columns or variants["SVLEN"].isnull().any():
        if not (variants.columns == "END").any() or variants["END"].isnull(
        ).any():
            variants["SVLEN"] = (
                variants["ALT"].str.len() - variants["REF"].str.len()).abs()
            print("REF ALT comp")
        else:
            print("use END")
            variants["SVLEN"] = variants["END"] - (variants["POS"] + 1)
    # convert to positive value
    variants.loc[:, "SVLEN"] = variants["SVLEN"].abs()
    if minlen is not None and maxlen is not None:
        variants = variants[(variants["SVLEN"] >= minlen)
                            & (variants["SVLEN"] < maxlen)]

    # only autosomes
    variants = variants[variants["CHROM"].str.match(r"(chr)?[0-9]+")]

    if constrain is not None:
        valid = (variants["MATCHING"] < 0) | (variants["MATCHING"].isin(
            constrain.index))
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
    if p == 0:
        return 1.0
    tp = np.count_nonzero(calls.is_tp)
    precision = tp / p
    return precision


def recall(calls, truth):
    p = calls.shape[0]
    if p == 0:
        return 0.0
    matches = calls.loc[calls.MATCHING.isin(truth.index), "MATCHING"]
    #tp = calls[calls.is_tp].MATCHING.unique().size
    tp = matches.unique().size
    t = truth.shape[0]
    recall = tp / t
    return recall


def get_colors(config):
    callers = [caller for caller in config["caller"] if caller != "varlociraptor"]
    palette = sns.color_palette("colorblind", n_colors=len(callers))
    palette = sns.color_palette("tab10", n_colors=len(callers))
    return {caller: c for caller, c in zip(callers, palette)}


def phred_scale(prob):
    return -10 * math.log10(prob)
