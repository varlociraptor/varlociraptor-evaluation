pindel_types = ["D", "SI"]#, "LI"]

rule pindel_config:
    input:
        get_bams
    output:
        "pindel/{run}.config.txt"
    run:
        with open(output[0], "w") as out:
            for f, t in zip(input, tissues):
                print(f, config["datasets"][config["runs"][wildcards.run]["dataset"]]["isize"]["mean"], t, file=out)


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
        extra=config["caller"]["pindel"]["params"]
    log:
        "logs/pindel/{run}.log"
    benchmark:
        "benchmarks/pindel/{run}.tsv"
    threads: 4
    wrapper:
        "0.19.1/bio/pindel/call"


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
        calls=expand("pindel/{{run}}.{vartype}.bcf", vartype=pindel_types),
        idx=expand("pindel/{{run}}.{vartype}.bcf.csi", vartype=pindel_types)
    output:
        "pindel/{run}.all.bcf"
    params:
        "-a -Ob"
    wrapper:
        "0.19.2/bio/bcftools/concat"
