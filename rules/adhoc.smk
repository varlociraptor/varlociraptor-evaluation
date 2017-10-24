rule adhoc_filter:
    input:
        vcf="{caller}/{run}.INDEL.bcf",
        bams=bams
    output:
        "adhoc-calls/{run}.bcf"
    conda:
        "envs/cyvcf2.yaml"
    script:
        "scripts/adhoc-calling.py"
