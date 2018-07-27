def get_prosic_input(ext):
    def get_prosic_input(wildcards):
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
    return get_prosic_input


rule prosic_call:
    input:
        calls=get_prosic_input(".bcf"),
        idx=get_prosic_input(".bcf.csi"),
        ref=get_ref,
        bams=get_bams,
        bais=get_bais
    output:
        temp("prosic-{caller}/{run}.{chrom}.bcf")
    params:
        caller=lambda wc: config["caller"]["prosic"].get(wc.caller, ""),
        chrom_prefix=lambda wc: config["ref"][config["runs"][wc.run]["ref"]].get("chrom_prefix", "") + wc.chrom,
        purity=lambda wc: config["runs"][wc.run]["purity"]
    log:
        "logs/prosic-{caller}/{run}.{chrom}.log"
    benchmark:
        "benchmarks/prosic-{caller}/{run}.{chrom}.tsv"
    # conda:
    #     "../envs/prosic.yaml"
    shell:
        "bcftools view {input.calls} {params.chrom_prefix} | "
        "prosic call-tumor-normal {input.bams} {input.ref} "
        "--purity {params.purity} "
        "{config[caller][prosic][params]} {params.caller} "
        "> {output} 2> {log}"


rule prosic_merge:
    input:
        expand("prosic-{{caller}}/{{run}}.{chrom}.bcf", chrom=CHROMOSOMES)
    output:
        "prosic-{caller}/{run}.all.bcf"
    wrapper:
        "0.19.1/bio/bcftools/concat"


rule prosic_control_fdr:
    input:
        "prosic-{caller}/{run}.all.bcf"
    output:
        "prosic-{caller}/{run}.gamma.{type}.{minlen}-{maxlen}.tsv"
    # conda:
    #     "../envs/prosic.yaml"
    shell:
        "prosic control-fdr --event SOMATIC --var {wildcards.type} "
        "--min-len {wildcards.minlen} --max-len {wildcards.maxlen} "
        "--method ev {input} > {output}"
