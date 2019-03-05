from itertools import product
from collections import namedtuple
from operator import itemgetter



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


rule match_varlociraptor_calls:
    input:
        calls="varlociraptor-{caller}/{run}.{vartype}.{minlen}-{maxlen}.{fdr}.bcf",
        truth=get_truth
    output:
        "matched-calls/varlociraptor-{caller}/{run}.{vartype}.{minlen}-{maxlen}.{fdr}.bcf"
    params:
        config["vcf-match-params"]
    #conda:
    #    "../envs/rbt.yaml"
    shell:
        "rbt vcf-match {params} {input.truth} < {input.calls} > {output}"


rule match_other_calls:
    input:
        calls="{mode}-{caller}/{run}.all.bcf",
        truth=get_truth
    output:
        "matched-calls/{mode}-{caller}/{run}.all.bcf"
    params:
        config["vcf-match-params"]
    #conda:
    #    "../envs/rbt.yaml"
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


def get_bcf_tags(wildcards):
    mode = wildcards.get("mode", "varlociraptor")
    caller = wildcards.caller
    if mode == "varlociraptor":
        caller = "varlociraptor"
    info_tags = list(config["caller"][caller].get("info", []))
    fmt_tags = list(config["caller"][caller].get("fmt", []))

    tags = ""
    if fmt_tags:
        tags += " --fmt {}".format(" ".join(fmt_tags))
    if info_tags:
        tags += " --info {}".format(" ".join(info_tags))
    return tags


def get_genotypes_param(wildcards):
    mode = wildcards.get("mode", "varlociraptor")
    caller = wildcards.caller
    if mode == "varlociraptor":
        caller = "varlociraptor"
    return "--genotypes" if config["caller"][caller].get("genotypes") else ""


rule varlociraptor_calls_to_tsv:
    input:
        "matched-calls/varlociraptor-{caller}/{run}.{vartype}.{minlen}-{maxlen}.{fdr}.bcf"
    output:
        "matched-calls/varlociraptor-{caller}/{run}.{vartype}.{minlen}-{maxlen}.{fdr}.tsv"
    params:
        tags=get_bcf_tags,
        gt=get_genotypes_param
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-to-txt {params.gt} {params.tags} --info MATCHING < {input} > {output}"


rule varlociraptor_all_calls_to_tsv:
    input:
        "varlociraptor-{caller}/{run}.all.bcf"
    output:
        "varlociraptor-{caller}/{run}.all.tsv"
    params:
        tags=get_bcf_tags,
        gt=get_genotypes_param
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-to-txt {params.gt} {params.tags} < {input} > {output}"


rule other_calls_to_tsv:
    input:
        "matched-calls/{mode}-{caller}/{run}.all.bcf"
    output:
        "matched-calls/{mode}-{caller}/{run}.all.tsv"
    params:
        tags=get_bcf_tags,
        gt=get_genotypes_param
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-to-txt {params.gt} {params.tags} --info MATCHING < {input} > {output}"


def get_tsv_calls(wildcards):
    pattern = "matched-calls/{mode}-{caller}/{run}.all.tsv"
    if wildcards.mode == "varlociraptor":
        pattern = "matched-calls/{mode}-{caller}/{run}.{vartype}.{minlen}-{maxlen}.{fdr}.tsv"
    return pattern.format(**wildcards)


rule obtain_tp_fp:
    input:
        calls=get_tsv_calls
    output:
        "annotated-calls/{mode}-{caller}/{run}.{vartype}.{minlen}-{maxlen}.{fdr}.tsv"
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/obtain-tp-fp.py"


def get_callers(mode):
    if mode == "varlociraptor":
        blacklist = config["caller"]["varlociraptor"]["blacklist"]
        callers = [caller for caller in config["caller"] if caller != "varlociraptor" and caller not in blacklist]
    elif mode == "default":
        callers = [caller for caller, p in config["caller"].items() if "score" in p and caller != "varlociraptor"]
    elif mode == "adhoc":
        callers = [caller for caller, p in config["caller"].items() if p.get("adhoc", False) and caller != "varlociraptor"]
    else:
        raise ValueError("Invalid mode: " + mode)
    return callers


def get_caller_runs(mode, runs):
    callers = get_callers(mode)
    return [(c, r) for c, r in product(callers, runs)]


def get_len_ranges(wildcards):
    return config["len-ranges"][wildcards.vartype]


def get_calls(mode, runs=None, fdr=[1.0]):
    def inner(wildcards):
        caller_runs = get_caller_runs(mode, [wildcards.run] if not runs else runs)
        len_ranges = get_len_ranges(wildcards)

        pattern = "annotated-calls/{mode}-{caller_run[0]}/{caller_run[1]}.{vartype}.{len_range[0]}-{len_range[1]}.{fdr}.tsv"
        return expand(pattern, mode=mode, caller_run=caller_runs,
                      vartype=wildcards.vartype, len_range=len_ranges, fdr=fdr)

    return inner


rule plot_precision_recall:
    input:
        varlociraptor_calls=get_calls("varlociraptor"),
        default_calls=get_calls("default"),
        adhoc_calls=get_calls("adhoc"),
        truth=lambda wc: "truth/{dataset}.annotated.tsv".format(**config["runs"][wc.run])
    output:
        "plots/precision-recall/{run}.{vartype}.svg"
    params:
        varlociraptor_callers=get_callers("varlociraptor"),
        default_callers=get_callers("default"),
        adhoc_callers=get_callers("adhoc"),
        len_ranges=get_len_ranges
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/plot-precision-recall.py"


rule plot_allelefreq_recall:
    input:
        varlociraptor_calls=get_calls("varlociraptor"),
        adhoc_calls=get_calls("adhoc"),
        truth=lambda wc: "truth/{dataset}.annotated.tsv".format(**config["runs"][wc.run])
    output:
        "plots/allelefreq-recall/{run}.{vartype}.svg"
    params:
        varlociraptor_callers=get_callers("varlociraptor"),
        adhoc_callers=get_callers("adhoc"),
        len_ranges=get_len_ranges
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/plot-allelefreq-recall.py"



rule plot_fdr:
    input:
        varlociraptor_calls=get_calls("varlociraptor", fdr=alphas),
    output:
        "plots/fdr-control/{run}.{vartype}.svg"
    params:
        callers=get_callers("varlociraptor"),
        purity=lambda wc: config["runs"][wc.run]["purity"],
        len_ranges=get_len_ranges,
        fdrs=alphas
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/plot-fdr-control.py"


rule plot_allelefreq:
    input:
        varlociraptor_calls=get_calls("varlociraptor"),
        truth=lambda wc: "truth/{dataset}.annotated.tsv".format(**config["runs"][wc.run])
    output:
        "plots/allelefreqs/{run}.{vartype}.svg"
    params:
        varlociraptor_callers=get_callers("varlociraptor"),
        len_ranges=get_len_ranges
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/plot-allelefreq-estimation.py"


rule plot_score_dist:
    input:
        varlociraptor_calls=get_calls("varlociraptor")
    output:
        "plots/score-dist/{run}.{vartype}.svg"
    params:
        varlociraptor_callers=get_callers("varlociraptor"),
        len_ranges=get_len_ranges
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/plot-score-dist.py"


rule plot_softclips:
    input:
        get_bams
    output:
        "plots/softclips/{run}.svg"
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/bam-stats.py"


rule concordance_match:
    input:
        lambda wc: expand("{mode}-{caller}/{run}.all.bcf" if wc.mode != "varlociraptor" else "varlociraptor-{caller}/{run}.adhoc.{threshold}.bcf",
                          run=config["plots"]["concordance"][wc.id], **wc),
    output:
        "concordance/{mode}-{caller}-{threshold}/{id}.{i}-vs-{j}.bcf"
    params:
        match=config["vcf-match-params"],
        bcfs=lambda wc, input: (input[int(wc.i)], input[int(wc.j)])
    #conda:
    #    "../envs/rbt.yaml"
    shell:
        "rbt vcf-match {params.match} {params.bcfs[1]} < {params.bcfs[0]} > {output}"


def get_concordance_combinations(id):
    assert len(config["plots"]["concordance"][id]) == 4
    return [(0, 1), (1, 2), (2, 0), (2, 3), (3, 0), (3, 1)]


rule aggregate_concordance:
    input:
        calls=lambda wc: expand("concordance/{mode}-{caller}-{threshold}/{id}.{i[0]}-vs-{i[1]}.tsv", i=get_concordance_combinations(wc.id), **wc),
        varlociraptor_calls=lambda wc: expand("varlociraptor-{caller}/{dataset}.all.tsv",
                                       dataset=[config["plots"]["concordance"][wc.id][c[0]] for c in get_concordance_combinations(wc.id)], **wc)
    output:
        "aggregated-concordance/{mode}-{caller}-{threshold}/{id}.{vartype}.tsv"
    params:
        dataset_combinations=lambda wc: list(get_concordance_combinations(wc.id)),
        datasets=lambda wc: config["plots"]["concordance"][wc.id]
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/aggregate-concordance.py"


rule concordance_to_tsv:
    input:
        "concordance/{mode}-{caller}-{threshold}/{prefix}.bcf"
    output:
        "concordance/{mode}-{caller}-{threshold}/{prefix}.tsv"
    params:
        tags=get_bcf_tags,
        gt=get_genotypes_param
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-to-txt {params.gt} {params.tags} --info MATCHING < {input} > {output}"



rule plot_concordance_upset:
    input:
        "aggregated-concordance/{mode}-{caller}-{threshold}/{id}.{vartype}.tsv"
    output:
        "plots/concordance/upset/{mode}-{caller}-{threshold}/{id}.{vartype}.concordance-upset.svg"
    params:
         datasets=lambda wc: config["plots"]["concordance"][wc.id]
    conda:
        "../envs/upset.yaml"
    script:
        "../scripts/plot-concordance-upset.R"


concordance_varlociraptor_calls = lambda threshold: expand("aggregated-concordance/varlociraptor-{caller}-{threshold}/{{id}}.{{vartype}}.tsv", caller=non_varlociraptor_callers, threshold=threshold)


rule plot_concordance:
    input:
        varlociraptor_calls_low=concordance_varlociraptor_calls(0.90),
        varlociraptor_calls_high=concordance_varlociraptor_calls(0.99),
        adhoc_calls=expand("aggregated-concordance/adhoc-{caller}-default/{{id}}.{{vartype}}.tsv", caller=non_varlociraptor_callers)
    output:
        "plots/concordance/{id}.{vartype}.concordance.svg"
    params:
        callers=non_varlociraptor_callers,
        len_ranges=get_len_ranges,
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/plot-concordance.py"


rule lancet_frequencies:
    input:
        "default-lancet/COLO_829-GSC.all.bcf"
    output:
        "plots/freqdist/lancet/COLO_829-GSC.svg"
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/plot-freqdist.py"
