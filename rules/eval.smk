rule annotate_truth:
    input:
        lambda wc: config["datasets"][wc.dataset]["truth"]
    output:
        "truth/{dataset}.annotated.bcf"
    conda:
        "../envs/cyvcf2.yaml"
    script:
        "../scripts/annotate-truth.py"


def get_truth(wildcards):
    if wildcards.mode == "prosic":
        run, _, purity = wildcards.run.rpartition("-")
    else:
        run = wildcards.run
    return "truth/{dataset}.annotated.bcf".format(**config["runs"][run])
        


rule match_variants:
    input:
        calls="{mode}-{caller}/{run}.all.bcf",
        truth=get_truth
    output:
        "matched-calls/{mode}-{caller}/{run}.all.bcf"
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-match {params} {input.truth} < {input.calls} > {output}"


rule truth_to_tsv:
    input:
        "truth/{dataset}.annotated.bcf"
    output:
        "truth/{dataset}.annotated.tsv"
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-to-txt --info SOMATIC SVLEN SVTYPE AF < {input} > {output}"


def get_info_tags(wildcards):
    tags = config["caller"][wildcards.caller].get("info", [])
    if wildcards.mode == "prosic":
        tags += config["caller"][wildcards.mode].get("info", [])
    return tags


rule calls_to_tsv:
    input:
        "matched-calls/{mode}-{caller}/{run}.all.bcf"
    output:
        "matched-calls/{mode}-{caller}/{run}.all.tsv"
    params:
        info=get_info_tags
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-to-txt --genotypes --info {params.info} MATCHING SVLEN SVTYPE < {input} > {output}"
