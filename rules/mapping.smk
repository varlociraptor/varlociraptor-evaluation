def get_source_bam(wildcards):
    ds = config["datasets"][wildcards.dataset]
    tissue = ds[wildcards.tissue]
    return tissue["bam"]


rule prepare_bam:
    input:
        get_source_bam
    output:
        temp("reads/{dataset}.{tissue}.namesorted.bam")
    params:
        "-n -m 2G"
    threads: 64
    wrapper:
        "0.22.0/bio/samtools/sort"


rule bam2fq:
    input:
        "reads/{dataset}.{tissue}.namesorted.bam"
    output:
        m1="reads/{dataset}.{tissue}.1.fastq.gz",
        m2="reads/{dataset}.{tissue}.2.fastq.gz",
        mixed=temp("reads/{dataset}.{tissue}.fastq.gz")
    conda:
        "../envs/tools.yaml"
    shell:
        "samtools bam2fq {input} -1 {output.m1} -2 {output.m2} -0 {output.mixed}"


rule bwa_index:
    input:
        lambda wc: config["ref"][wc.ref]["fasta"]
    output:
        "index/{ref}/genome.bwt"
    params:
        prefix=lambda wc, output: os.path.splitext(output[0])[0]
    conda:
        "../envs/bwa.yaml"
    shell:
        "bwa index -p {params.prefix} {input}"


bwa_params = r'-Y -R "@RG\\tID:{tissue}\\tSM:{tissue}"'


def get_reads(wildcards):
    ds = config["datasets"][wildcards.dataset][wildcards.tissue]
    if "fq1" in ds:
        return [ds["fq1"], ds["fq2"]]
    else:
        return expand("reads/{{dataset}}.{{tissue}}.{mate}.fastq.gz", mate=[1, 2])


rule bwa:
    input:
        index="index/{ref}/genome.bwt",
        sample=get_reads
    output:
        temp("mapped-bwa/{dataset}.{tissue}.{ref}.bam")
    log:
        "logs/bwa/{dataset}.{tissue}.{ref}.log"
    benchmark:
        "benchmarks/bwa/{dataset}.{tissue}.{ref}.tsv"
    params:
        index=lambda wc, input: os.path.splitext(input.index)[0],
        extra=bwa_params
    conda:
        "../envs/bwa.yaml"
    threads: 8
    shell:
        "(resources/bwa mem -t {threads} {params.extra} {params.index} {input.sample} | "
        "samtools view -Sb - > {output}) 2> {log}"


rule samtools_sort:
    input:
        "mapped-{mapper}/{dataset}.{tissue}.{ref}.bam"
    output:
        temp("mapped-{mapper}/{dataset}.{tissue}.{ref}.sorted.pre.bam")
    params:
        "-m 2G"
    threads: 64
    wrapper:
        "0.22.0/bio/samtools/sort"


rule mark_duplicates:
    input:
        "mapped-{mapper}/{dataset}.{tissue}.{ref}.sorted.pre.bam"
    output:
        bam=protected("mapped-{mapper}/{dataset}.{tissue}.{ref}.sorted.bam"),
        metrics="mapped-{mapper}/{dataset}.{tissue}.{ref}.markdup.metrics.txt"
    params:
        "-Xmx6g VALIDATION_STRINGENCY=LENIENT -Djava.io.tmpdir=tmp"
    log:
        "logs/picard/dedup/{mapper}/{dataset}.{tissue}.{ref}.log"
    wrapper:
        "0.22.0/bio/picard/markduplicates"


rule samtools_index:
    input:
        "mapped-{mapper}/{dataset}.{tissue}.{ref}.sorted.bam"
    output:
        "mapped-{mapper}/{dataset}.{tissue}.{ref}.sorted.bam.bai"
    wrapper:
        "0.22.0/bio/samtools/index"

