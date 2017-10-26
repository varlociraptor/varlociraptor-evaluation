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
        "-a"
    wrapper:
        "0.17.4/bio/bcftools/concat"


rule delly_samples:
    output:
        "resources/{dataset}.delly-samples.txt"
    run:
        with open(output[0], "w") as out:
            ds = config["datasets"][wildcards.dataset]
            print(ds["tumor"]["name"], "tumor", sep="\t", file=out)
            print(ds["normal"]["name"], "control", sep="\t", file=out)


ruleorder: delly_adhoc > adhoc_filter


rule delly_adhoc:
    input:
        bcf="delly/{run}.all.bcf",
        samples=lambda wc: "resources/{dataset}.delly-samples.txt".format(
            dataset=config["runs"][wc.run]["dataset"])
    output:
        "adhoc-calls/delly/{run}.bcf"
    conda:
        "../envs/delly.yaml"
    shell:
        "delly filter -m 0 -r 1.0 --samples {input.samples} "
        "-o {output} {input.bcf}"
