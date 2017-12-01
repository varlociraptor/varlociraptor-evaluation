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


def get_targets(run):
    non_prosic = [caller for caller in config["caller"] if caller != "prosic" and caller != "pindel"]
    t = expand("adhoc-{caller}/{run}.all.bcf",
               caller=non_prosic,
               run=run)
    t.extend(expand("lancet/{run}.all.bcf", run=run))
    t.extend(expand("prosic-{caller}/{run}-{purity}.all.bcf",
                    run=run, caller=non_prosic,
                    purity=config["runs"][run]["purity"]))
    return t


wildcard_constraints:
    chrom="|".join(CHROMOSOMES),
    caller="|".join(config["caller"]),
    tissue="tumor|normal",
    ref="|".join(config["ref"]),
    run="|".join(config["runs"]),
    purity="[01]\.[0-9]+"


rule all:
    input:
        [get_targets(run) for run in config["runs"] if run != "test"]


rule test:
    input:
        get_targets("test")


include: "rules/mapping.smk"
include: "rules/delly.smk"
include: "rules/pindel.smk"
include: "rules/lancet.smk"
include: "rules/prosic.smk"
include: "rules/adhoc.smk"


rule index_bcf:
    input:
        "{prefix}.bcf"
    output:
        "{prefix}.bcf.csi"
    conda:
        "envs/bcftools.yaml"
    shell:
        "bcftools index {input}"
