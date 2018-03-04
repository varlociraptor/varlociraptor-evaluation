rule strelka:
    input:
        ref=get_ref,
        samples=get_bams,
        bais=get_bais,
        manta="manta/{run}/results/variants/candidateSmallIndels.vcf.gz"
    output:
        "strelka/{run}/results/variants/somatic.snvs.vcf.gz",
        "strelka/{run}/results/variants/somatic.indels.vcf.gz"
    params:
        dir=lambda w, output: os.path.dirname(os.path.dirname(os.path.dirname(output[0])))
    log:
        "logs/strelka/{run}.log"
    benchmark:
        "benchmarks/strelka/{run}.tsv"
    conda:
        "../envs/strelka.yaml"
    threads: 16
    shell:
        "rm -rf {params.dir}; "
        "(configureStrelkaSomaticWorkflow.py --tumorBam {input.samples[0]} --normalBam {input.samples[1]} "
        "--referenceFasta {input.ref} --runDir {params.dir} --indelCandidates {input.manta}; "
        "{params.dir}/runWorkflow.py -m local -j {threads}) > {log} 2>&1"


rule strelka_default:
    input:
        "strelka/{run}/results/variants/somatic.indels.vcf.gz",
        "strelka/{run}/results/variants/somatic.snvs.vcf.gz"
    output:
        "default-strelka/{run}.all.bcf"
    params:
        "-a -Ob"
    wrapper:
        "0.22.0/bio/bcftools/concat"


ruleorder: strelka_adhoc > adhoc_filter


rule strelka_adhoc:
    input:
        "default-strelka/{run}.all.bcf"
    output:
        "adhoc-strelka/{run}.all.bcf"
    params:
        "-i INFO/SOMATIC -f PASS -Ob"
    wrapper:
        "0.22.0/bio/bcftools/view"

