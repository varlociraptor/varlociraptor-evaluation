def get_whole_chrom_region(wildcards, input):
    chrom = "chr" + wildcards.chrom
    return chrom


rule lancet:
    input:
        ref=get_ref,
        bams=get_bams,
        bais=get_bais
    output:
        temp("lancet/{run}/chr{chrom}.vcf")
    params:
        region=get_whole_chrom_region,
        extra=config["caller"]["lancet"]["params"]
    log:
        "logs/lancet/{run}.chr{chrom}.log"
    benchmark:
        "benchmarks/lancet/{run}.chr{chrom}.tsv"
    wildcard_constraints:
        chrom="[^.]+"
    threads: 4
    shell:
       "resources/lancet --tumor {input.bams[0]} --normal {input.bams[1]} "
       "--ref {input.ref} --reg {params.region} --num-threads {threads} "
       "{params.extra} > {output} 2> {log}"


rule fix_lancet:
    input:
        header="resources/lancet_header.txt",
        vcf="lancet/{run}/chr{chrom}.vcf"
    output:
        temp("lancet/{run}/chr{chrom}.fixed.vcf")
    conda:
        "../envs/tools.yaml"
    shell:
        "sed -r 's/MS\=[0-9]+[ACGT]+/MS/g' {input.vcf} | bcftools annotate -o {output} -h {input.header} -"


rule merge_lancet:
    input:
        expand("lancet/{{run}}/chr{chrom}.fixed.vcf", chrom=CHROMOSOMES)
    output:
        "default-lancet/{run}.all.bcf"
    conda:
        "../envs/tools.yaml"
    shell:
        "bcftools concat -Ob {input} > {output}"
