def get_prosic_input(ext):
    def get_prosic_input(wildcards):
        return "{prefix}{caller}/{run}.all{ext}".format(
            prefix="" if "score" not in config["caller"][wildcards.caller] else "default-",
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
        temp("prosic-{caller}/{run}-{purity}.{chrom}.bcf")
    params:
        isize=lambda wc: config["datasets"][config["runs"][wc.run]["dataset"]]["isize"],
        caller=lambda wc: config["caller"]["prosic"].get(wc.caller, ""),
        chrom_prefix=lambda wc: config["ref"][config["runs"][wc.run]["ref"]].get("chrom_prefix", "") + wc.chrom
    log:
        "logs/prosic-{caller}/{run}-{purity}.{chrom}.log"
    benchmark:
        "benchmarks/prosic-{caller}/{run}-{purity}.{chrom}.tsv"
    # conda:
    #     "../envs/prosic.yaml"
    shell:
        "bcftools view {input.calls} {params.chrom_prefix} | "
        "prosic call-tumor-normal {input.bams} {input.ref} "
        "--isize-mean {params.isize[mean]} --isize-sd {params.isize[sd]} "
        "--purity {wildcards.purity} "
        "{config[caller][prosic][params]} {params.caller} "
        "> {output} 2> {log}"


rule prosic_merge:
    input:
        expand("prosic-{{caller}}/{{run}}-{{purity}}.{chrom}.bcf", chrom=CHROMOSOMES)
    output:
        "prosic-{caller}/{run}-{purity}.all.bcf"
    wrapper:
        "0.19.1/bio/bcftools/concat"


rule prosic_control_fdr:
    input:
        "prosic-{caller}/{run}-{purity}.all.bcf"
    output:
        "prosic-{caller}/{run}-{purity}.gamma.{type}.{minlen}-{maxlen}.tsv"
    # conda:
    #     "../envs/prosic.yaml"
    shell:
        "prosic control-fdr --event SOMATIC --var {wildcards.type} "
        "--min-len {wildcards.minlen} --max-len {wildcards.maxlen} "
        "--method ev {input} > {output}"
