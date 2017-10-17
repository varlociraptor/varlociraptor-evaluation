

rule pindel_config:
    input:
        bams
    output:
        "pindel/{dataset}.config.txt"
    run:
        with open(output[0], "w") as out:
            for f, t in zip(input, tissues):
                print(f, 312, t)


rule pindel:
    input:
        ref="index/hg19.fa",
        # samples to call
        samples=bams,
        # bam configuration file, see http://gmt.genome.wustl.edu/packages/pindel/quick-start.html
        config="pindel/{dataset}.config.txt"
    output:
        expand("pindel/{dataset}_{type}", type=pindel_types)
    params:
        # prefix must be consistent with output files
        prefix="pindel/{dataset}",
        extra=""  # optional parameters (except -i, -f, -o)
    log:
        "logs/pindel.log"
    threads: 4
    wrapper:
        "0.17.4/bio/pindel/call"



rule pindel2bcf:
    input:
        ref="index/hg19.fa",
        pindel="pindel/{dataset}_{type}",
        header="resources/contigs.vcf"
    output:
        "pindel/{dataset}.{type}.bcf"
    params:
        refname=config["ref"]["name"],  # mandatory, see pindel manual
        refdate=config["ref"]["date"],  # mandatory, see pindel manual
    log:
        "logs/pindel/{sample}/{sample}.{contig}.pindel2bcf.{type}.log"
    conda:
        "../envs/pindel2bcf.yaml"
    shell:
        "(pindel2vcf -p {input.pindel} -r {input.ref} -R {params.refname} "
        "-d {params.refdate} -v {output}.vcf && "
        "bcftools annotate -Ou -h {input.header} {output}.vcf | "
        "bcftools view -Ob -l9 -o {output} - && rm {output}.vcf) "
        "> {log} 2>&1"
