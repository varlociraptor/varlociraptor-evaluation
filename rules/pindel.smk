pindel_types = ["D", "SI", "LI"]

rule pindel_config:
    input:
        get_bams
    output:
        "pindel/{run}.config.txt"
    run:
        with open(output[0], "w") as out:
            for f, t in zip(input, tissues):
                print(f, 312, t)


rule pindel:
    input:
        ref=get_ref,
        # samples to call
        samples=get_bams,
        bais=get_bais,
        # bam configuration file, see http://gmt.genome.wustl.edu/packages/pindel/quick-start.html
        config="pindel/{run}.config.txt"
    output:
        expand("pindel/{{run}}_{type}", type=pindel_types)
    params:
        # prefix must be consistent with output files
        prefix="pindel/{run}",
        extra=config["caller"]["delly"]["params"]
    log:
        "logs/pindel/{run}.log"
    threads: 4
    wrapper:
        "0.17.4/bio/pindel/call"


rule pindel2bcf:
    input:
        ref=get_ref,
        pindel="pindel/{run}_{type}",
        header="resources/contigs.vcf"
    output:
        "pindel/{run}.{type}.bcf"
    params:
        refname=lambda wc: config["runs"][wc.run]["ref"],
        refdate=lambda wc: config["ref"][config["runs"][wc.run]["ref"]]["date"]
    log:
        "logs/pindel/{run}.{type}.log"
    conda:
        "../envs/pindel.yaml"
    shell:
        "(pindel2vcf -p {input.pindel} -r {input.ref} -R {params.refname} "
        "-d {params.refdate} -v {output}.vcf && "
        "bcftools annotate -Ou -h {input.header} {output}.vcf | "
        "bcftools view -Ob -l9 -o {output} - && rm {output}.vcf) "
        "> {log} 2>&1"


rule pindel_concat:
    input:
        expand("pindel/{{run}}.{vartype}.bcf", vartype=pindel_types)
    output:
        "pindel/{run}.all.bcf"
    params:
        "-a"
    wrapper:
        "0.17.4/bio/bcftools/concat"
