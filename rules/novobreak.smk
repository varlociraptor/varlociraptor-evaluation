rule novobreak:
    input:
        ref=get_ref,
        bams=get_bams,
        bais=get_bais
    output:
        workdir=directory("novobreak/{run}")
    log:
        "logs/novobreak/{run}.log"
    benchmark:
        "benchmarks/novobreak/{run}.tsv"
    conda:
        "../envs/novobreak.yaml"
    threads: 16
    shell:
        "run_novobreak {input.ref} {input.bams[0]} {input.bams[1]} {threads} {output.workdir}"


rule novobreak_adhoc:
    input:
        "novobreak/{run}"
    output:
        "adhoc-novobreak/{run}.all.bcf"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools view -Ob {input}/novoBreak.pass.flt.vcf > {output}"


rule novobreak_all:
    input:
        "novobreak/{run}"
    output:
        "novobreak/{run}.all.bcf"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools concat -Ob {input}/split/*.sp.vcf > {output}"


ruleorder: neusomatic_adhoc > adhoc_filter

