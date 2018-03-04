rule scalpel_discovery:
    input:
        ref=get_ref,
        samples=get_bams,
        bais=get_bais
    output:
        "scalpel/{run}/main/database.db",
        "scalpel/{run}/twopass/database.db"
    params:
        dir=lambda wc, output: os.path.dirname(output[0])
    log:
        "logs/scalpel-discovery/{run}.log"
    benchmark:
        "benchmarks/scalpel-discovery/{run}.tsv"
    conda:
        "../envs/scalpel.yaml"
    shell:
        "rm -rf {params.dir}; "
        "scalpel-discovery --somatic --normal {input[1]} --tumor {input[0]} "
        "--ref {input.ref} --dir {params.dir} --two-pass 2> {log} 2>&1"


rule scalpel_export:
    input:
        ref=get_ref,
        db="scalpel/{run}/{mode}/database.db"
    output:
        "scalpel/{run}/{mode}/variants.all.vcf"
    params:
        dir=lambda wc, output: os.path.dirname(output[0])
    conda:
        "../envs/scalpel.yaml"
    shell:
        "scalpel-export --somatic --db {input.db} --output-format vcf --variant-type all "
        "--ref {input.ref} --dir {params.dir} 2> {log} 2>&1"
