import pandas as pd

from snakemake.remote import FTP
from snakemake.remote import EGA

ega = EGA.RemoteProvider()

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


def get_stats(wildcards):
    run = config["runs"][wildcards.run]
    return expand("stats-{mapper}/{dataset}.{tissue}.{ref}.stats.txt",
                  dataset=run["dataset"],
                  tissue=tissues,
                  ref=run["ref"],
                  mapper=run["mapper"])


def get_ref(wildcards):
    ref = config["runs"][wildcards.run]["ref"]
    return config["ref"][ref]["fasta"]


def get_caller_params(caller):
    def inner(wildcards):
        return "{} {}".format(config["caller"][caller].get("params", ""),
                              config["runs"][wildcards.run].get("params", dict()).get(caller, ""))
    return inner


wildcard_constraints:
    chrom="|".join(CHROMOSOMES),
    caller="|".join(config["caller"]),
    tissue="tumor|normal",
    ref="|".join(config["ref"]),
    run="|".join(config["runs"]),
    purity="[01]\.[0-9]+",
    mode="prosic|adhoc|default"


target_concordance = expand("plots/concordance/{id}.{vartype}.{lenrange[0]}-{lenrange[1]}.concordance.svg", id=config["plots"]["concordance"], vartype=["INS", "DEL"], lenrange=config["len-ranges"])


rule all:
    input:
        expand("plots/{plt}/{run}.{vartype}.{lenrange[0]}-{lenrange[1]}.svg",
               plt=["precision-recall", "fdr-control"],
               lenrange=config["len-ranges"],
               vartype=["INS", "DEL"],
               run=config["plots"]["known-truth"]),
        expand("plots/allelefreqs/{run}.{vartype}.{lenrange[0]}-{lenrange[1]}.svg",
               lenrange=config["len-ranges"],
               vartype=["INS", "DEL"],
               run=config["plots"]["known-truth"]),
        target_concordance

rule all_concordance:
    input:
        target_concordance


include: "rules/mapping.smk"
include: "rules/delly.smk"
#include: "rules/pindel.smk"
include: "rules/lancet.smk"
include: "rules/manta.smk"
include: "rules/strelka.smk"
include: "rules/bpi.smk"
include: "rules/prosic.smk"
include: "rules/adhoc.smk"
include: "rules/eval.smk"
include: "rules/stats.smk"


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
        "ranked-{type,[ft]}ps/prosic-{caller}/{run}-0.75.{vartype}.{minlen}-{maxlen}.tsv"
    params:
        type=lambda wc: False if wc.type == "f" else True
    run:
        calls = pd.read_table(input[0])
        calls = calls[calls.is_tp == params.type]
        calls.sort_values("PROB_SOMATIC", ascending=not params.type, inplace=True)
        calls[["CHROM", "POS", "PROB_SOMATIC", "is_tp"]].to_csv(output[0], sep="\t")


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
    wildcard_constraints:
        caller=".+"
    params:
        region=testcase_region
    conda:
        "envs/bcftools.yaml"
    shell:
        "bcftools view {input.bcf} {wildcards.varpos} > {output.vcf}; "
        "samtools view -b {input.bams[0]} {params.region} > {output.tumor}; "
        "samtools view -b {input.bams[1]} {params.region} > {output.normal}; "
        "samtools index {output.tumor}; samtools index {output.normal}"
