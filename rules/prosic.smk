rule prosic_call:
    input:
        calls="{caller}/{run}.all.bcf",
        idx="{caller}/{run}.all.bcf.csi",
        ref=get_ref,
        bams=get_bams,
        bais=get_bais
    output:
        temp("prosic-{caller}/{run}-{purity}.{chrom}.bcf")
    params:
        isize=lambda wc: config["datasets"][config["runs"][wc.run]["dataset"]]["isize"],
        caller=lambda wc: config["caller"]["prosic"].get(wc.caller, ""),
        chrom_prefix=lambda wc: config["refs"][config["runs"][wc.run]["ref"]].get("chrom_prefix", "")
    log:
        "logs/prosic-{caller}/{run}-{purity}.log"
    benchmark:
        "benchmarks/prosic-{caller}/{run}-{purity}.tsv"
    # conda:
    #     "../envs/prosic.yaml"
    shell:
        "bcftools view {input.calls} {wildcards.chrom} | "
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
        "prosic-{caller}/{run}-{purity}.gamma.{type}.tsv"
    # conda:
    #     "../envs/prosic.yaml"
    shell:
        "prosic control-fdr --event SOMATIC --var {wildcards.type} "
        "--method ev {input} > {output}"
