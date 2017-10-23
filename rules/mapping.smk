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
        lambda wc: config["datasets"][wc.dataset][wc.tissue]
    output:
        "reads/{dataset}.{tissue}.1.fastq",
        "reads/{dataset}.{tissue}.2.fastq"
    shell:
        "samtools bam2fq {input} {output}"


rule bowtie2_index:
    input:
        ftp.remote(config["ref"]["index"], keep_local=True, static=True)
    output:
        "index"
    shell:
        "mkdir -p index; tar -C {output} -xf {input}"


rule qtip:
    input:
        index="index",
        m1="reads/{dataset}.{tissue}.1.fastq",
        m2="reads/{dataset}.{tissue}.2.fastq"
    output:
        "mapped/{dataset}.{tissue}.bam"
    conda:
        "envs/qtip.yaml"
    log:
        "logs/qtip/{dataset}.{tissue}.log"
    threads: 8
    shell:
        "(qtip --bt2-exe 'bowtie2 -p {threads}' "
        "--m1 {input.m1} --m2 {input.m2} --index {input.index} | "
        "samtools view -Sb - > {output} && samtools index {output}) 2> {log}"


rule bowtie2:
    input:
        sample=expand("reads/{{dataset}}.{{tissue}}.{mate}.fastq", mate=[1, 2])
    output:
        "mapped/{dataset}.{tissue}.bam"
    log:
        "logs/bowtie2/{dataset}.{tissue}.log"
    params:
        index="index/hg38",  # prefix of reference genome index (built with bowtie2-build)
        extra=""  # optional parameters
    threads: 8
    wrapper:
        "0.17.4/bio/bowtie2/align"
