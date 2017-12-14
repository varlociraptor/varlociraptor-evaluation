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
    if wildcards.mode == "prosic":
        run, _, purity = wildcards.run.rpartition("-")
    else:
        run = wildcards.run
    return "truth/{dataset}.annotated.{ext}".format(ext=ext, **config["runs"][run])


rule match_variants:
    input:
        calls="{mode}-{caller}/{run}.all.bcf",
        truth=get_truth
    output:
        "matched-calls/{mode}-{caller}/{run}.all.bcf"
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
        "rbt vcf-to-txt --genotypes --info SOMATIC SVLEN SVTYPE AF < {input} > {output}"


def get_info_tags(wildcards):
    tags = list(config["caller"][wildcards.caller].get("info", []))
    if wildcards.mode == "prosic":
        tags += config["caller"][wildcards.mode].get("info", [])
    return tags


rule calls_to_tsv:
    input:
        "matched-calls/{mode}-{caller}/{run}.all.bcf"
    output:
        "matched-calls/{mode}-{caller}/{run}.all.tsv"
    params:
        info=get_info_tags
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-to-txt --genotypes --info {params.info} MATCHING SVLEN SVTYPE < {input} > {output}"


rule obtain_tp_fp:
    input:
        truth=partial(get_truth, ext="tsv"),
        calls="matched-calls/{mode}-{caller}/{run}.all.tsv"
    output:
        "matched-calls/{mode}-{caller}/{run}.{vartype}.{minlen}-{maxlen}.tsv"
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/obtain-tp-fp.py"


def get_callers(mode):
    if mode == "prosic":
        callers = [caller for caller in config["caller"] if caller != "prosic"]
    elif mode == "default":
        callers = [caller for caller, p in config["caller"].items() if "score" in p and caller != "prosic"]
    elif mode == "adhoc":
        callers = [caller for caller, p in config["caller"].items() if "score" not in p]
    else:
        raise ValueError("Invalid mode: " + mode)
    return callers


def get_calls(mode):
    callers = get_callers(mode)

    def inner(wildcards):
        if mode == "prosic":
            purity = config["runs"][wildcards.run]["purity"]
            sep = "-"
        else:
            purity = ""
            sep = ""
        return expand("matched-calls/{mode}-{caller}/{run}{sep}{purity}.{vartype}.{minlen}-{maxlen}.tsv", mode=mode, caller=callers, purity=purity, sep=sep, **wildcards)

    return inner


rule plot_precision_recall:
    input:
        prosic_calls=get_calls("prosic"),
        default_calls=get_calls("default"),
        adhoc_calls=get_calls("adhoc"),
        truth=lambda wc: "truth/{dataset}.annotated.tsv".format(**config["runs"][wc.run])
    output:
        "plots/precision-recall/{run}.{vartype}.{minlen}-{maxlen}.svg"
    params:
        prosic_callers=get_callers("prosic"),
        default_callers=get_callers("default"),
        adhoc_callers=get_callers("adhoc"),
        purity=lambda wc: config["runs"][wc.run]["purity"]
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/plot-precision-recall.py"
