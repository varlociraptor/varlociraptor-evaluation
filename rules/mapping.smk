rule get_bam:
    input:
        config["gdc-manifest"]
    output:
        temp(expand("gdc-data/{bam}", bam=gdc_manifest.filename))
    conda:
        "envs/gdc-client.yaml"
    shell:
        "gdc-client download --debug -m {input}"


rule prepare_bam:
    input:
        lambda wc: config["datasets"][wc.dataset][wc.tissue]["bam"]
    output:
        temp("reads/{dataset}.{tissue}.namesorted.bam")
    params:
        "-n"
    threads: 8 
    wrapper:
        "0.18.0/bio/samtools/sort"


rule bam2fq:
    input:
        "reads/{dataset}.{tissue}.namesorted.bam"
    output:
        m1="reads/{dataset}.{tissue}.1.fastq",
        m2="reads/{dataset}.{tissue}.2.fastq"
        #mixed=temp("reads/{dataset}.{tissue}.fastq")
    conda:
        "../envs/tools.yaml"
    shell:
        "samtools bam2fq {input} -1 {output.m1} -2 {output.m2} "
        #"cat {output.mixed} | grep '^@.*/1$' -A 3 --no-group-separator > {output.m1} && "
        #"cat {output.mixed} | grep '^@.*/2$' -A 3 --no-group-separator > {output.m2}"


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
        ref=lambda wc: config["ref"][wc.ref]["fasta"],
        m1="reads/{dataset}.{tissue}.1.fastq",
        m2="reads/{dataset}.{tissue}.2.fastq"
    output:
        "mapped-qtip/{dataset}.{tissue}.{ref}.bam"
    params:
        index="index/{ref}/genome",
        tmp="mapped-qtip"
    conda:
        "../envs/qtip.yaml"
    log:
        "logs/qtip/{dataset}.{tissue}.{ref}.log"
    threads: 12
    shell:
        "(qtip --bt2-exe 'bowtie2 --local -p {threads}' --temp-directory {params.tmp} "
        "--m1 {input.m1} --m2 {input.m2} --index {params.index} --ref {input.ref} | "
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
