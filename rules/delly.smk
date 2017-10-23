rule delly:
    input:
        ref="index/hg19.fa",
        samples=bams,
    output:
        "delly/{run}.{type,(DEL|DUP|INV|TRA|INS)}.bcf"
    params:
        vartype="{type}", # variant type to call (can be wildcard, hardcoded string or function)
        extra=""  # optional parameters for delly (except -t, -g)
    log:
        "logs/delly/{run}.{type}.log"
    threads: 2
    wrapper:
        "0.17.4/bio/delly"


rule delly_concat:
    input:
        expand("delly/{{run}}.{type}.bcf", vartype=["DEL", "INS"])
    output:
        "delly/{run}.INDEL.bcf"
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


rule delly_adhoc:
    input:
        bcf="delly/{run}.INDEL.bcf",
        samples=lambda wc: "resources/{dataset}.delly-samples.txt".format(
            dataset=config["runs"][wc.run]["dataset"])
    output:
        "adhoc-calls/delly/{run}.bcf"
    conda:
        "envs/delly.yaml"
    shell:
        "delly filter -m 0 -r 1.0 --samples {input.samples} -o {output} {input.bcf}"
