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
        lambda wc: ftp.remote(config["ref"][wc.ref]["index"], keep_local=True, static=True)
    output:
        "index/{ref}/genome.fasta"
    params:
        folder="index/{ref}"
    shell:
        "tar -C {params.folder} -xf {input}"


rule qtip:
    input:
        index="index/{ref}",
        m1="reads/{dataset}.{tissue}.1.fastq",
        m2="reads/{dataset}.{tissue}.2.fastq"
    output:
        "mapped/{dataset}.{tissue}.{ref}.bam"
    conda:
        "envs/qtip.yaml"
    log:
        "logs/qtip/{dataset}.{tissue}.{ref}.log"
    threads: 8
    shell:
        "(qtip --bt2-exe 'bowtie2 -p {threads}' "
        "--m1 {input.m1} --m2 {input.m2} --index {input.index} | "
        "samtools view -Sb - > {output} && samtools index {output}) 2> {log}"


rule bowtie2:
    input:
        index="index/{ref}",
        sample=expand("reads/{{dataset}}.{{tissue}}.{mate}.fastq", mate=[1, 2])
    output:
        "mapped/{dataset}.{tissue}.{ref}.bam"
    log:
        "logs/bowtie2/{dataset}.{tissue}.{ref}.log"
    params:
        index="index/{ref}",  # prefix of reference genome index
        extra=""  # optional parameters
    threads: 8
    wrapper:
        "0.17.4/bio/bowtie2/align"
