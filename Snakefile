import pandas as pd

from snakemake.remote import FTP

ftp = FTP.RemoteProvider()


configfile: "config.yaml"


gdc_manifest = pd.read_table(config["gdc-manifest"])
units = pd.DataFrame({
    "sample": gdc_manifest.filename.str.slice(5, 33),
    "path": expand("gdc-data/{bam}", bam=gdc_manifest.filename)})
units.index = gdc_manifest.id


rule all:
    input:
        expand("mapped/{unit}.bam", unit=units.index)


rule get_bam:
    input:
        config["gdc-manifest"]
    output:
        temp(expand("gdc-data/{bam}", bam=gdc_manifest.filename))
    conda:
        "envs/gdc-client.yaml"
    shell:
        "gdc-client download --debug -m {input}"


rule bam2fq:
    input:
        lambda wc: config["datasets"][wc.dataset][wc.tissue]
    output:
        "reads/{dataset}.{tissue}.1.fastq",
        "reads/{dataset}.{tissue}.2.fastq"
    shell:
        "samtools bam2fq {input} {output}"


rule bowtie2_index:
    input:
        ftp.remote(config["ref"]["index"], keep_local=True, static=True)
    output:
        "index"
    shell:
        "mkdir -p index; tar -C {output} -xf {input}"


rule qtip:
    input:
        index="index",
        m1="reads/{dataset}.{tissue}.1.fastq",
        m2="reads/{dataset}.{tissue}.2.fastq"
    output:
        "mapped/{dataset}.{tissue}.bam"
    conda:
        "envs/qtip.yaml"
    threads: 8
    shell:
        "qtip --bt2-exe 'bowtie2 -p {threads}' "
        "--m1 {input.m1} --m2 {input.m2} --index {input.index} | "
        "samtools view -Sb - > {output} && samtools index {output}"


# TODO calling and concordance analysis

############# callers ####################

tissues = ["tumor", "normal"]
bams = expand("mapped/{{dataset}}.{tissue}.bam", tissue=tissues)


rule delly:
    input:
        ref="index/hg19.fa",
        samples=bams,
    output:
        "delly/{dataset}.{type,(DEL|DUP|INV|TRA|INS)}.bcf"
    params:
        vartype="{type}", # variant type to call (can be wildcard, hardcoded string or function)
        extra=""  # optional parameters for delly (except -t, -g)
    log:
        "logs/delly/{dataset}.{type}.log"
    threads: 2
    wrapper:
        "0.17.4/bio/delly"


rule delly_concat:
    input:
        expand("delly/{{dataset}}.{type}.bcf", vartype=["DEL", "INS"])
    output:
        "delly/{dataset}.INDEL.bcf"
    params:
        "-a"
    wrapper:
        "0.17.4/bio/bcftools/concat"


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





########### adhoc calling #############


rule delly_adhoc:
    input:
        bcf="external-fixed/simulated-delly.bcf",
        samples="resources/delly-samples.txt"
    output:
        "results/adhoc-calls/simulated-delly.vcf"
    params:
        temp="results/adhoc-calls/simulated-delly.bcf"
    conda:
        "envs/delly.yaml"
    shell:
        "delly filter -m 0 -r 1.0 --samples {input.samples} -o {params.temp} {input.bcf} "
        "&& bcftools view {params.temp} > {output}"


rule 
