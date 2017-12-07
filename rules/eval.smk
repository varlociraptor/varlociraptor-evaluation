rule annotate_truth:
    input:
        lambda wc: config["datasets"][wc.dataset]["truth"]
    output:
        "truth/{dataset}.annotated.bcf"
    conda:
        "../envs/cyvcf2.yaml"
    script:
        "../scripts/annotate-truth.py"


rule match_variants:
    input:
        calls="{mode}-{caller}/{run}.all.bcf",
        truth="truth/{dataset}.annotated.bcf"
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
    shell:
        "rbt vcf-to-txt --info SOMATIC SVLEN SVTYPE AF < {input} > {output}"


rule calls_to_tsv:
    input:
        "{mode}-{caller}/{run}.all.bcf"
    output:
        "{mode}-{caller}/{run}.all.tsv"
    params:
        info=lambda wc: " ".join(config["caller"][wc.caller].get("info", []))
    shell:
        "rbt vcf-to-txt --genotypes --info {params.info} MATCHING SVLEN SVTYPE < {input} > {output}"
