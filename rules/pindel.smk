pindel_types = ["D", "SI", "LI"]

rule pindel_config:
    input:
        bams
    output:
        "pindel/{run}.config.txt"
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
        config="pindel/{run}.config.txt"
    output:
        expand("pindel/{run}_{type}", type=pindel_types)
    params:
        # prefix must be consistent with output files
        prefix="pindel/{run}",
        extra=""  # optional parameters (except -i, -f, -o)
    log:
        "logs/pindel/{run}.log"
    threads: 4
    wrapper:
        "0.17.4/bio/pindel/call"


rule pindel2bcf:
    input:
        ref="index/hg19.fa",
        pindel="pindel/{run}_{type}",
        header="resources/contigs.vcf"
    output:
        "pindel/{run}.{type}.bcf"
    params:
        refname=config["ref"]["name"],  # mandatory, see pindel manual
        refdate=config["ref"]["date"],  # mandatory, see pindel manual
    log:
        "logs/pindel/{run}.{type}.log"
    conda:
        "../envs/pindel2bcf.yaml"
    shell:
        "(pindel2vcf -p {input.pindel} -r {input.ref} -R {params.refname} "
        "-d {params.refdate} -v {output}.vcf && "
        "bcftools annotate -Ou -h {input.header} {output}.vcf | "
        "bcftools view -Ob -l9 -o {output} - && rm {output}.vcf) "
        "> {log} 2>&1"


rule pindel_concat:
    input:
        expand("pindel/{{run}}.{type}.bcf", vartype=["DEL", "INS"])
    output:
        "pindel/{run}.all.bcf"
    params:
        "-a"
    wrapper:
        "0.17.4/bio/bcftools/concat"
