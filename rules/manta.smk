rule manta:
    input:
        ref=get_ref,
        samples=get_bams,
        bais=get_bais
    output:
        "manta/{run}/results/variants/candidateSV.vcf.gz",
        "manta/{run}/results/variants/somaticSV.vcf.gz",
        "manta/{run}/results/variants/candidateSmallIndels.vcf.gz"
    params:
        dir=lambda w, output: os.path.dirname(os.path.dirname(os.path.dirname(output[0])))
    log:
        "logs/manta/{run}.log"
    benchmark:
        "benchmarks/manta/{run}.tsv"
    conda:
        "../envs/manta.yaml"
    threads: 8
    shell:
        "configManta.py --tumorBam {input.samples[0]} --normalBam {input.samples[1]} "
        "--referenceFasta {input.ref} --runDir {params.dir}; "
        "{params.dir}/runWorkflow.py -m local -j {threads}"


rule manta_default:
    input:
        "manta/{run}/results/variants/somaticSV.vcf.gz"
    output:
        "default-manta/{run}.all.bcf"
    params:
        "-Ob"
    wrapper:
        "0.22.0/bio/bcftools/view"


ruleorder: manta_adhoc > adhoc_filter


rule manta_adhoc:
    input:
        "default-manta/{run}.all.bcf"
    output:
        "adhoc-manta/{run}.all.bcf"
    params:
        "-f PASS -i INFO/SOMATIC -Ob"
    wrapper:
        "0.22.0/bio/bcftools/view"

