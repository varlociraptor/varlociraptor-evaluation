rule prepare_bam:
    input:
        lambda wc: config["datasets"][wc.dataset][wc.tissue]["bam"]
    output:
        temp("reads/{dataset}.{tissue}.namesorted.bam")
    params:
        "-n"
    threads: 8 
    wrapper:
        "0.22.0/bio/samtools/sort"


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


rule bwa_index:
    input:
        lambda wc: config["ref"][wc.ref]["fasta"]
    output:
        "index/{ref}/genome.bwt"
    params:
        prefix=lambda wc, output: os.path.splitext(output[0])[0]
    conda:
        "../envs/qtip.yaml"
    shell:
        "bwa index -p {params.prefix} {input}"


bwa_params = r'-Y -R "@RG\\tID:{tissue}\\tSM:{tissue}"'


rule qtip:
    input:
        index="index/{ref}/genome.bwt",
        ref=lambda wc: config["ref"][wc.ref]["fasta"],
        m1="reads/{dataset}.{tissue}.1.fastq",
        m2="reads/{dataset}.{tissue}.2.fastq"
    output:
        temp("mapped-qtip/{dataset}.{tissue}.{ref}.bam")
    params:
        index=lambda wc, input: os.path.splitext(input.index)[0],
        out="mapped-qtip/{dataset}.{tissue}.{ref}",
        tmp="mapped-qtip",
        bwa=bwa_params
    conda:
        "../envs/qtip.yaml"
    log:
        "logs/qtip/{dataset}.{tissue}.{ref}.log"
    benchmark:
        "benchmarks/qtip/{dataset}.{tissue}.{ref}.tsv"
    threads: 8
    resources:
        disk=1  # constrain number of qtip instances
    shell:
        "rm -rf {params.out}; "
        "(resources/qtip --bwa-exe 'resources/bwa mem {params.bwa} -t {threads}' "
        "--output-directory {params.out} --temp-directory {params.tmp} --input-model-size 30000 "
        "--aligner bwa-mem --m1 {input.m1} --m2 {input.m2} --index {params.index} --ref {input.ref}; "
        "samtools view -Sb {params.out}/final.sam > {output}; "
        "rm -r {params.out}) 2> {log}"


rule bwa:
    input:
        index="index/{ref}/genome.bwt",
        sample=expand("reads/{{dataset}}.{{tissue}}.{mate}.fastq", mate=[1, 2])
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
        "../envs/qtip.yaml"
    threads: 8
    shell:
        "(resources/bwa mem -t {threads} {params.extra} {params.index} {input.sample} | "
        "samtools view -Sb - > {output}) 2> {log}"


rule samtools_sort:
    input:
        "mapped-{mapper}/{dataset}.{tissue}.{ref}.bam"
    output:
        temp("mapped-{mapper}/{dataset}.{tissue}.{ref}.sorted.pre.bam")
    threads: 8
    wrapper:
        "0.22.0/bio/samtools/sort"


rule mark_duplicates:
    input:
        "mapped-{mapper}/{dataset}.{tissue}.{ref}.sorted.pre.bam"
    output:
        protected("mapped-{mapper}/{dataset}.{tissue}.{ref}.sorted.bam")
    params:
        "-Xmx2g VALIDATION_STRINGENCY=LENIENT"
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


