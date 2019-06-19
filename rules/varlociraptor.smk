import math


def get_varlociraptor_input(ext):
    def get_varlociraptor_input(wildcards):
        if wildcards.caller == "truth":
            return "truth/{dataset}.annotated.bcf".format(
                dataset=config["runs"][wildcards.run]["dataset"])

        conf = config["caller"][wildcards.caller]
        prefix = ""
        if "score" in conf and not conf.get("useraw"):
            prefix = "default-"
        return "{prefix}{caller}/{run}.all{ext}".format(
            prefix=prefix,
            ext=ext, **wildcards)
    return get_varlociraptor_input


rule varlociraptor_call:
    input:
        calls=get_varlociraptor_input(".bcf"),
        idx=get_varlociraptor_input(".bcf.csi"),
        ref=get_ref,
        bams=get_bams,
        bais=get_bais
    output:
        temp("varlociraptor-{caller}/{run}.{chrom}.bcf")
    params:
        caller=lambda wc: config["caller"]["varlociraptor"].get(wc.caller, ""),
        chrom_prefix=lambda wc: config["ref"][config["runs"][wc.run]["ref"]].get("chrom_prefix", "") + wc.chrom,
        purity=lambda wc: config["runs"][wc.run]["purity"]
    log:
        "logs/varlociraptor-{caller}/{run}.{chrom}.log"
    benchmark:
        "benchmarks/varlociraptor-{caller}/{run}.{chrom}.tsv"
    # conda:
    #     "../envs/varlociraptor.yaml"
    shell:
        "bcftools view -Ou {input.calls} {params.chrom_prefix} | "
        "varlociraptor call variants {input.ref} "
        "{config[caller][varlociraptor][params]} {params.caller} "
        "tumor-normal {input.bams} "
        "--purity {params.purity} "
        "> {output} 2> {log}"


rule varlociraptor_merge:
    input:
        expand("varlociraptor-{{caller}}/{{run}}.{chrom}.bcf", chrom=CHROMOSOMES)
    output:
        "varlociraptor-{caller}/{run}.all.bcf"
    params:
        "-Ob"
    wrapper:
        "0.19.1/bio/bcftools/concat"


rule varlociraptor_filter_by_odds:
    input:
        "varlociraptor-{caller}/{run}.all.bcf"
    output:
        "varlociraptor-{caller}/{run}.oddsfiltered.bcf"
    shell:
        "varlociraptor filter-calls posterior-odds positive --events SOMATIC_TUMOR < {input} > {output}"


rule varlociraptor_control_fdr:
    input:
        "varlociraptor-{caller}/{run}.all.bcf"
    output:
        "varlociraptor-{caller}/{run}.{type}.{minlen}-{maxlen}.{fdr}.bcf"
    # conda:
    #     "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor filter-calls control-fdr {input} --events SOMATIC_TUMOR --var {wildcards.type} "
        "--minlen {wildcards.minlen} --maxlen {wildcards.maxlen} "
        "--fdr {wildcards.fdr} > {output}"


rule adhoc_varlociraptor:
    input:
        "varlociraptor-{caller}/{run}.all.bcf"
    output:
        "varlociraptor-{caller}/{run}.adhoc.{threshold}.bcf"
    params:
        lambda wc: "-i 'INFO/PROB_SOMATIC_TUMOR<={}' -Ob".format(-10 * math.log10(float(wc.threshold)))
    wildcard_constraints:
        threshold="0.[0-9]+"
    wrapper:
        "0.22.0/bio/bcftools/view"
