rule delly:
    input:
        ref=get_ref,
        samples=get_bams,
        bais=get_bais
    output:
        "delly/{run}.{type,(DEL|DUP|INV|TRA|INS)}.bcf"
    params:
        vartype="{type}", # variant type to call
        extra=config["caller"]["delly"]["params"]
    log:
        "logs/delly/{run}.{type}.log"
    benchmark:
        "benchmarks/delly/{run}.{type}.tsv"
    threads: 2
    wrapper:
        "0.17.4/bio/delly"


rule delly_concat:
    input:
        expand("delly/{{run}}.{vartype}.bcf", vartype=["DEL", "INS"])
    output:
        "delly/{run}.all.bcf"
    params:
        "-a -Ob"
    wrapper:
        "0.17.4/bio/bcftools/concat"


ruleorder: delly_adhoc > adhoc_filter


rule delly_adhoc:
    input:
        csi="delly/{run}.all.bcf.csi",
        bcf="delly/{run}.all.bcf",
        samples="resources/delly-samples.txt"
    output:
        "adhoc-calls/delly/{run}.bcf"
    conda:
        "../envs/delly.yaml"
    shell:
        "delly filter -m 0 -r 1.0 --samples {input.samples} "
        "-o {output} {input.bcf}"
