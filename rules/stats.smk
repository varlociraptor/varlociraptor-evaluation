rule insert_size_dist:
    """Calculate insert size distribution."""
    input:
        lambda wc: config["datasets"][wc.dataset][wc.tissue]["bam"]
    output:
        txt="tables/insert-size/{dataset}.{tissue}.txt",
        pdf="plots/insert-size/{dataset}.{tissue}.pdf"
    log:
        "logs/picard/insert-size/{dataset}.{tissue}.log"
    params:
        "VALIDATION_STRINGENCY=LENIENT "
        "METRIC_ACCUMULATION_LEVEL=null "
        "METRIC_ACCUMULATION_LEVEL=SAMPLE "
    wrapper:
        "0.22.0/bio/picard/collectinsertsizemetrics"
