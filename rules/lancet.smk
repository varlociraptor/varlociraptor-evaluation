rule build_lancet:
    output:
        "resources/lancet"
    conda:
        "../envs/lancet.yaml"
    params:
        commit="c3d08f58a681279539ed55178ae42a3625905cd5"
    shell:
        """
        cd resources
        curl -L https://github.com/nygenome/lancet/archive/{params.commit}.tar.gz | tar -xz
        cd lancet-{params.commit}
        sed -i 's/-Wl,-rpath,$(ABS_BAMTOOLS_DIR)\/lib\///' Makefile
        make clean
        prefix=$(dirname $(dirname $(which $CXX)))
        make CXX=$CXX INCLUDES="-I$prefix/include/ -I$prefix/include/bamtools -L$prefix/lib"
        cp lancet ..
        """

def get_whole_chrom_region(wildcards, input):
    if config["runs"][wildcards.run]["ref"] == "hg18":
        chrom = "chr" + wildcards.chrom
    else:
        chrom = wildcards.chrom
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
        extra=get_caller_params("lancet")
    log:
        "logs/lancet/{run}.chr{chrom}.log"
    benchmark:
        "benchmarks/lancet/{run}.chr{chrom}.tsv"
    wildcard_constraints:
        chrom="[^.]+"
    conda:
        "../envs/lancet.yaml"
    threads: 4
    shell:
       "LD_LIBRARY_PATH=$CONDA_PREFIX/lib "
       "resources/lancet --tumor {input.bams[0]} --normal {input.bams[1]} "
       "--ref {input.ref} --reg {params.region} --num-threads {threads} "
       "{params.extra} > {output} 2> {log}"


rule fix_lancet:
    input:
        header=lambda w: "resources/lancet.{}.header.txt".format(config["runs"][w.run]["ref"]),
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


ruleorder: lancet_adhoc > adhoc_filter


rule lancet_adhoc:
    input:
        "default-lancet/{run}.all.bcf"
    output:
        "adhoc-lancet/{run}.all.bcf"
    params:
        "-Ob -f PASS -i INFO/SOMATIC"
    wrapper:
        "0.19.3/bio/bcftools/view"
