rule samtools_stats:
    """Calculate insert size distribution."""
    input:
        "mapped-{mapper}/{dataset}.{tissue}.{ref}.sorted.bam"
    output:
        "stats-{mapper}/{dataset}.{tissue}.{ref}.stats.txt"
    log:
        "logs/samtools/stats/{mapper}.{dataset}.{tissue}.{ref}.log"
    wrapper:
        "0.22.0/bio/samtools/stats"
