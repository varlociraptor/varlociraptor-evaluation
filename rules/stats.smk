rule samtools_stats:
    """Calculate insert size distribution."""
    input:
        bam="mapped-{mapper}/{dataset}.{tissue}.{ref}.sorted.bam"
    output:
        "stats-{mapper}/{dataset}.{tissue}.{ref}.stats.txt"
    log:
        "logs/samtools/stats/{mapper}.{dataset}.{tissue}.{ref}.log"
    params:
        region="2:1000-20000000"
    wrapper:
        "0.22.0/bio/samtools/stats"


rule depth:
    input:
        "mapped-{mapper}/{dataset}.{tissue}.{ref}.sorted.bam"
    output:
        "stats-{mapper}/{dataset}.{tissue}.{ref}.depth.per-base.tsv.gz"
    params:
        prefix=lambda w, output: output[:-len(".per-base.bed.gz")]
    conda:
        "../envs/bcftools.yaml"
    shell:
        "samtools depth {input} | gzip > {output}"
