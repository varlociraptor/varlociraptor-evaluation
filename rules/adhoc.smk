rule adhoc_filter:
    input:
        vcf="{caller}/{run}.all.bcf",
        bams=get_bams
    output:
        "adhoc-{caller}/{run}.all.bcf"
    conda:
        "../envs/cyvcf2.yaml"
    script:
        "scripts/adhoc-calling.py"
