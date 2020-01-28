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
non_varlociraptor_callers = [caller for caller in config["caller"] if caller != "varlociraptor" and caller != "pindel"]
alphas = [0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
vartypes = ["INS", "DEL"]

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
    mode="varlociraptor|adhoc|default",
    minlen="[0-9]+",
    maxlen="[0-9]+"




rule all:
    input:
        expand("plots/{plt}/{run}.{vartype}.{ext}",
               plt=["fdr-control", "allelefreqs", "score-dist", "allelefreq-recall", "allelefreq-scatter"],
               vartype=vartypes,
               run=config["plots"]["known-truth-full-afs"],
               ext=["pdf"]),
        expand("plots/precision-recall/{run}.{vartype}.{zoom}.{ext}",
               vartype=vartypes,
               run=config["plots"]["known-truth"],
               zoom=["zoom", "nozoom"],
               ext=["pdf"]),
        expand("plots/concordance/{run}.{vartype}.concordance.{ext}",
               run=config["plots"]["concordance"],
               vartype=vartypes,
               ext=["pdf"])


include: "rules/mapping.smk"
include: "rules/delly.smk"
#include: "rules/pindel.smk"
include: "rules/lancet.smk"
include: "rules/manta.smk"
include: "rules/strelka.smk"
include: "rules/neusomatic.smk"
include: "rules/novobreak.smk"
include: "rules/bpi.smk"
include: "rules/varlociraptor.smk"
include: "rules/adhoc.smk"
include: "rules/eval.smk"
include: "rules/stats.smk"
include: "rules/tools.smk"


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
        "matched-calls/varlociraptor-{caller}/{run}.{vartype}.{minlen}-{maxlen}.1.0.tsv"
    output:
        "ranked-{type,[ft]}ps/varlociraptor-{caller}/{run}-0.75.{vartype}.{minlen}-{maxlen}.tsv"
    params:
        type=lambda wc: False if wc.type == "f" else True
    run:
        calls = pd.read_table(input[0], header=[0, 1])
        calls = calls["VARIANT"]
        calls = calls[(calls.MATCHING >= 0 if params.type else calls.MATCHING < 0)]
        calls.sort_values("PROB_SOMATIC_TUMOR", ascending=not params.type, inplace=True)
        calls[["CHROM", "POS", "PROB_SOMATIC_TUMOR", "MATCHING"]].to_csv(output[0], sep="\t")


def testcase_region(wildcards):
    chrom, pos = wildcards.varpos.split(":")
    pos = int(pos)
    return "{}:{}-{}".format(chrom, pos - 1000, pos + 1000)


rule testcase:
    input:
        bcf=get_varlociraptor_input(".bcf"),
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
