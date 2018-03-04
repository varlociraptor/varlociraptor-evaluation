rule bpi:
    input:
        manta="manta/{run}/results/variants/somaticSV.vcf.gz",
        samples=get_bams,
        bais=get_bais
    output:
        "bpi/{run}.all.vcf",
    log:
        "logs/bpi/{run}.log"
    benchmark:
        "benchmarks/bpi/{run}.tsv"
    conda:
        "../envs/bpi.yaml"
    shell:
        "break-point-inspector -vcf {input.manta} -ref {input.samples[1]} -tumor {input.samples[0]} "
        "-output_vcf {output}"


ruleorder: bpi_adhoc > adhoc_filter


rule bpi_adhoc:
    input:
        "bpi/{run}.all.vcf"
    output:
        "bpi/{run}.all.bcf"
    params:
        "-i INFO/SOMATIC -f PASS -Ob"
    wrapper:
        "0.22.0/bio/bcftools/view"
