rule prosic_call:
    input:
        calls="{caller}/{run}.all.bcf"
        ref=get_ref,
        bams=get_bams,
        bais=get_bais
    output:
        "prosic-{caller}/{run}.{chrom}.bcf"
    log:
        "logs/prosic-{caller}/{run}.log"
    benchmark:
        "benchmarks/prosic-{caller}/{run}.tsv"
    shell:
        "bcftools view {input.calls} {wildcards.chrom} | "
        "prosic call-tumor-normal {input.bams} {input.ref} "
        "{config[prosic][params]} > {output} 2> {log}"


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
        "prosic-{caller}/{run}.gamma.{type}.tsv"
    shell:
        "prosic control-fdr --event SOMATIC --var {wildcards.type} "
        "--method ev {input} > {output}"
