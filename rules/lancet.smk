def get_whole_chrom_region(wildcards, input):
    chrom = "chr" + wildcards.chrom
    idx = pd.read_table(input.ref + ".fai", index_col=0, names=["len", "offset", "linebases", "linewidth"])
    return "{}:1-{}".format(chrom, idx.loc[chrom][0])


rule lancet:
    input:
        ref=get_ref,
        bams=get_bams,
        bais=get_bais
    output:
        "lancet/{run}/chr{chrom}.vcf"
    params:
        region=get_whole_chrom_region
    log:
        "logs/lancet/{run}.chr{chrom}.log"
    wildcard_constraints:
        chrom="[^.]+"
    threads: 24
    shell:
       "lancet --tumor {input.bams[0]} --normal {input.bams[1]} "
       "--ref {input.ref} --reg {params.region} --num-threads {threads} "
       "> {output} 2> {log}"


rule fix_lancet:
    input:
        header="resources/lancet_header.txt",
        vcf="lancet/{run}/chr{chrom}.vcf"
    output:
        "lancet/{run}/chr{chrom}.fixed.vcf"
    conda:
        "envs/tools.yaml"
    shell:
        "sed -r 's/MS\=[0-9]+[ACGT]+/MS/g' {input.vcf} | bcftools annotate -o {output} -h {input.header} -"


rule merge_lancet:
    input:
        expand("lancet/{{run}}/chr{chrom}.fixed.vcf", chrom=CHROMOSOMES)
    output:
        "lancet/{run}/all.bcf"
    conda:
        "envs/tools.yaml"
    shell:
        "bcftools concat -Ob {input} > {output}"
