import pandas as pd

from snakemake.remote import FTP

ftp = FTP.RemoteProvider()


configfile: "config.yaml"


gdc_manifest = pd.read_table(config["gdc-manifest"])
units = pd.DataFrame({
    "sample": gdc_manifest.filename.str.slice(5, 33),
    "path": expand("gdc-data/{bam}", bam=gdc_manifest.filename)})
units.index = gdc_manifest.id
CHROMOSOMES = list(map(str, range(1,23))) + ["M", "X", "Y"]
tissues = ["tumor", "normal"]
non_prosic_callers = [caller for caller in config["caller"] if caller != "prosic" and caller != "pindel"]


def get_bams(wildcards):
    run = config["runs"][wildcards.run]
    return expand("mapped-{mapper}/{dataset}.{tissue}.{ref}.sorted.bam",
                  dataset=run["dataset"],
                  tissue=tissues,
                  ref=run["ref"],
                  mapper=run["mapper"])


def get_bais(wildcards):
    return ["{}.bai".format(f) for f in get_bams(wildcards)]


def get_ref(wildcards):
    ref = config["runs"][wildcards.run]["ref"]
    return config["ref"][ref]["fasta"]


wildcard_constraints:
    chrom="|".join(CHROMOSOMES),
    caller="|".join(config["caller"]),
    tissue="tumor|normal",
    ref="|".join(config["ref"]),
    #run="|".join(config["runs"]),
    purity="[01]\.[0-9]+",
    mode="prosic|adhoc|default"


rule all:
    input:
        expand("plots/{plt}/{run}.{vartype}.{lenrange[0]}-{lenrange[1]}.svg",
               plt=["precision-recall", "fdr-control"],
               lenrange=config["len-ranges"],
               vartype=["INS", "DEL"],
               run=config["runs"]),
        expand("plots/allelefreqs/{run[0]}-{run[1]}.{vartype}.{lenrange[0]}-{lenrange[1]}.svg",
               lenrange=config["len-ranges"],
               vartype=["INS", "DEL"],
               run=[(name, purity) for name, p in config["runs"].items() for purity in p["purity"]])


include: "rules/mapping.smk"
include: "rules/delly.smk"
#include: "rules/pindel.smk"
include: "rules/lancet.smk"
include: "rules/prosic.smk"
include: "rules/adhoc.smk"
include: "rules/eval.smk"


rule index_bcf:
    input:
        "{prefix}.bcf"
    output:
        "{prefix}.bcf.csi"
    conda:
        "envs/bcftools.yaml"
    shell:
        "bcftools index {input}"


rule rank_fps:
    input:
        "matched-calls/prosic-{caller}/{run}-0.75.{vartype}.{minlen}-{maxlen}.tsv"
    output:
        "ranked-fps/prosic-{caller}/{run}-0.75.{vartype}.{minlen}-{maxlen}.tsv"
    params:
        cols=lambda wc: "9,2,3,17" if wc.caller == "delly" else "2,3,8,16"
    shell:
        "set +o pipefail; cut -f {params.cols} {input} | grep False | sort -n -k3 > {output}"


def testcase_region(wildcards):
    chrom, pos = wildcards.varpos.split(":")
    pos = int(pos)
    return "{}:{}-{}".format(chrom, pos - 1000, pos + 1000)


rule testcase:
    input:
        bcf=get_prosic_input(".bcf"),
        bams=get_bams,
        bais=get_bais
    output:
        vcf="testcase/{run}/{caller}-{varpos}/candidates.vcf",
        tumor="testcase/{run}/{caller}-{varpos}/tumor.bam",
        normal="testcase/{run}/{caller}-{varpos}/normal.bam"
    params:
        region=testcase_region
    conda:
        "envs/bcftools.yaml"
    shell:
        "bcftools index -f {input.bcf}; "
        "bcftools view {input.bcf} {wildcards.varpos} > {output.vcf}; "
        "samtools view -b {input.bams[0]} {params.region} > {output.tumor}; "
        "samtools view -b {input.bams[1]} {params.region} > {output.normal}"
