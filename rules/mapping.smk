rule get_bam:
    input:
        config["gdc-manifest"]
    output:
        temp(expand("gdc-data/{bam}", bam=gdc_manifest.filename))
    conda:
        "envs/gdc-client.yaml"
    shell:
        "gdc-client download --debug -m {input}"


rule bam2fq:
    input:
        lambda wc: config["datasets"][wc.dataset][wc.tissue]["bam"]
    output:
        m1="reads/{dataset}.{tissue}.1.fastq",
        m2="reads/{dataset}.{tissue}.2.fastq"
    conda:
        "../envs/tools.yaml"
    shell:
        "samtools bam2fq {input} -1 {output.m1} -2 {output.m2}"


rule bowtie2_index:
    input:
        lambda wc: config["ref"][wc.ref]["fasta"]
    output:
        "index/{ref}/genome.1.bt2"
    params:
        prefix="index/{ref}/genome"
    conda:
        "../envs/qtip.yaml"
    shell:
        "bowtie2-build {input} {params.prefix}"


rule qtip:
    input:
        index="index/{ref}/genome.1.bt2",
        m1="reads/{dataset}.{tissue}.1.fastq",
        m2="reads/{dataset}.{tissue}.2.fastq"
    output:
        "mapped-qtip/{dataset}.{tissue}.{ref}.bam"
    params:
        index="index/{ref}/genome"
    conda:
        "../envs/qtip.yaml"
    log:
        "logs/qtip/{dataset}.{tissue}.{ref}.log"
    threads: 8
    shell:
        "(qtip --bt2-exe 'bowtie2 -p {threads}' "
        "--m1 {input.m1} --m2 {input.m2} --index {params.index} | "
        "samtools view -Sb - > {output}) 2> {log}"


rule bowtie2:
    input:
        index="index/{ref}/genome.1.bt2",
        sample=expand("reads/{{dataset}}.{{tissue}}.{mate}.fastq", mate=[1, 2])
    output:
        "mapped-bowtie2/{dataset}.{tissue}.{ref}.bam"
    log:
        "logs/bowtie2/{dataset}.{tissue}.{ref}.log"
    params:
        index="index/{ref}/genome",
        extra=""  # optional parameters
    threads: 8
    wrapper:
        "0.18.0/bio/bowtie2/align"


rule samtools_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    wrapper:
        "0.18.0/bio/samtools/index"
