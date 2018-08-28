from itertools import product
from collections import namedtuple


rule annotate_truth:
    input:
        lambda wc: config["datasets"][wc.dataset]["truth"]
    output:
        "truth/{dataset}.annotated.bcf"
    conda:
        "../envs/cyvcf2.yaml"
    script:
        "../scripts/annotate-truth.py"


def get_truth(wildcards, ext="bcf"):
    return "truth/{dataset}.annotated.{ext}".format(ext=ext, **config["runs"][wildcards.run])


rule match_variants:
    input:
        calls="{mode}-{caller}/{run}.{vartype}.{minlen}-{maxlen}.{fdr}.bcf",
        truth=get_truth
    output:
        "matched-calls/{mode}-{caller}/{run}.{vartype}.{minlen}-{maxlen}.{fdr}.bcf"
    params:
        config["vcf-match-params"]
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-match {params} {input.truth} < {input.calls} > {output}"


rule truth_to_tsv:
    input:
        "truth/{dataset}.annotated.bcf"
    output:
        "truth/{dataset}.annotated.tsv"
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-to-txt --genotypes --info SOMATIC SVLEN SVTYPE TAF NAF < {input} > {output}"


def get_info_tags(wildcards):
    tags = list(config["caller"][wildcards.caller].get("info", []))
    if wildcards.mode == "prosic":
        tags += config["caller"][wildcards.mode].get("info", [])
    return tags


rule calls_to_tsv:
    input:
        "matched-calls/{mode}-{caller}/{run}.{vartype}.{minlen}-{maxlen}.{fdr}.bcf"
    output:
        "matched-calls/{mode}-{caller}/{run}.{vartype}.{minlen}-{maxlen}.{fdr}.tsv"
    params:
        info=get_info_tags,
        gt=lambda wc: "--genotypes" if config["caller"][wc.caller].get("genotypes") else ""
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-to-txt {params.gt} --info {params.info} MATCHING < {input} > {output}"


rule obtain_tp_fp:
    input:
        calls="matched-calls/{mode}-{caller}/{run}.{vartype}.{minlen}-{maxlen}.{fdr}.tsv"
    output:
        "annotated-calls/{mode}-{caller}/{run}.{vartype}.{minlen}-{maxlen}.{fdr}.tsv"
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/obtain-tp-fp.py"


def get_callers(mode):
    if mode == "prosic":
        blacklist = config["caller"]["prosic"]["blacklist"]
        callers = [caller for caller in config["caller"] if caller != "prosic" and caller not in blacklist]
    elif mode == "default":
        callers = [caller for caller, p in config["caller"].items() if "score" in p and caller != "prosic"]
    elif mode == "adhoc":
        callers = [caller for caller, p in config["caller"].items() if p.get("adhoc", False) and caller != "prosic"]
    else:
        raise ValueError("Invalid mode: " + mode)
    return callers


def get_caller_runs(mode, runs):
    callers = get_callers(mode)
    return [(c, r) for c, r in product(callers, runs)]


def get_calls(mode, runs=None, fdr=[1.0], len_range=config["len-ranges"]):
    def inner(wildcards):
        caller_runs = get_caller_runs(mode, [wildcards.run] if not runs else runs)

        pattern = "annotated-calls/{mode}-{caller_run[0]}/{caller_run[1]}.{vartype}.{len_range[0]}-{len_range[1]}.{fdr}.tsv"
        return expand(pattern, mode=mode, caller_run=caller_runs,
                      vartype=wildcards.vartype, len_range=len_range, fdr=fdr)

    return inner


rule plot_precision_recall:
    input:
        prosic_calls=get_calls("prosic"),
        default_calls=get_calls("default"),
        adhoc_calls=get_calls("adhoc"),
        truth=lambda wc: "truth/{dataset}.annotated.tsv".format(**config["runs"][wc.run])
    output:
        "plots/precision-recall/{run}.{vartype}.svg"
    params:
        prosic_callers=get_callers("prosic"),
        default_callers=get_callers("default"),
        adhoc_callers=get_callers("adhoc"),
        len_ranges=config["len-ranges"]
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/plot-precision-recall.py"


rule plot_fdr:
    input:
        prosic_calls=get_calls("prosic", fdr=alphas),
    output:
        "plots/fdr-control/{run}.{vartype}.svg"
    params:
        props=list(product(get_callers("prosic"), alphas)),
        purity=lambda wc: config["runs"][wc.run]["purity"]
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/plot-fdr-control.py"


rule plot_allelefreq:
    input:
        prosic_calls=get_calls("prosic"),
        truth=lambda wc: "truth/{dataset}.annotated.tsv".format(**config["runs"][wc.run])
    output:
        "plots/allelefreqs/{run}.{vartype}.svg"
    params:
        prosic_callers=get_callers("prosic")
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/plot-allelefreq-estimation.py"


rule plot_softclips:
    input:
        get_bams
    output:
        "plots/softclips/{run}.svg"
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/bam-stats.py"


rule fig_fdr:
    input:
        expand("plots/fdr-control/{{run}}.{{vartype}}.{lenrange[0]}-{lenrange[1]}.svg",
               lenrange=config["len-ranges"])
    output:
        "figs/{run}.{vartype}.fdr.svg"
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/fig-fdr.py"


def get_concordance_calls(mode, files=True):
    def inner(wildcards):
        runs = config["plots"]["concordance"][wildcards.id]
        if files:
            callers = get_callers(mode)
            return expand("concordance/{mode}-{caller}/{id}.{i}.tsv",
                          mode=mode,
                          caller=callers,
                          id=wildcards.id,
                          i=[0, 1])
        else:
            return get_caller_runs(mode, runs)
    return inner


def get_depths(wildcards):
    """Returns depth files for the two runs defined by the given concordance id."""
    runs = [config["runs"][r] for r in config["plots"]["concordance"][wildcards.id]]
    return expand("stats-{run[mapper]}/{run[dataset]}.{tissue}.{run[ref]}.depth.per-base.tsv.gz",
                  run=runs,
                  tissue=tissues)


def get_concordance_bams(wildcards):
    """Returns bam files for the two runs defined by the given concordance id."""
    runs = [config["runs"][r] for r in config["plots"]["concordance"][wildcards.id]]
    return expand("mapped-{run[mapper]}/{run[dataset]}.{tissue}.{run[ref]}.sorted.bam",
                  run=runs,
                  tissue=tissues)


rule concordance_match:
    input:
        lambda wc: expand("{mode}-{caller}/{run}.all.bcf",
                          run=config["plots"]["concordance"][wc.id], **wc),
        exons="concordance/{id}.covered-exons.bed"
    output:
        "concordance/{mode}-{caller}/{id}.{i}.bcf"
    params:
        match=config["vcf-match-params"],
        bcfs=lambda wc, input: (input[int(wc.i)], input[1 - int(wc.i)])
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-match {params.match} {params.bcfs[0]} < {params.bcfs[1]} | "
        "bcftools view | "
        "bedtools intersect -header -b {input.exons} -a - -wa -f 1.0 | "
        "bcftools view -Ob > {output}"


rule concordance_to_tsv:
    input:
        "concordance/{mode}-{caller}/{id}.{i}.bcf"
    output:
        "concordance/{mode}-{caller}/{id}.{i}.tsv"
    params:
        info=get_info_tags,
        gt=lambda wc: "--genotypes" if config["caller"][wc.caller].get("genotypes") else ""
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-to-txt {params.gt} --info {params.info} MATCHING < {input} > {output}"


rule covered_exons:
    input:
        bams=get_concordance_bams,
        exons=lambda wildcards: config["ref"][config["runs"][config["plots"]["concordance"][wildcards.id][0]]["ref"]]["exons"]
    output:
        bed="concordance/{id}.covered-exons.bed"
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/exoncov.py"


rule plot_concordance:
    input:
        prosic_calls=get_concordance_calls("prosic"),
        default_calls=get_concordance_calls("default"),
        adhoc_calls=get_concordance_calls("adhoc")
    output:
        "plots/concordance/{id}.{vartype}.{minlen}-{maxlen}.concordance.svg"
    params:
        prosic_runs=get_concordance_calls("prosic", files=False),
        default_runs=get_concordance_calls("default", files=False),
        adhoc_runs=get_concordance_calls("adhoc", files=False),
        len_ranges=config["len-ranges"],
        runs=lambda wc: config["plots"]["concordance"][wc.id]
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/plot-concordance.py"
