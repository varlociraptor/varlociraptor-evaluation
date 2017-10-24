import pandas as pd

from snakemake.remote import FTP

ftp = FTP.RemoteProvider()


configfile: "config.yaml"


gdc_manifest = pd.read_table(config["gdc-manifest"])
units = pd.DataFrame({
    "sample": gdc_manifest.filename.str.slice(5, 33),
    "path": expand("gdc-data/{bam}", bam=gdc_manifest.filename)})
units.index = gdc_manifest.id


def get_bams(wildcards):
    tissues = ["tumor", "normal"]
    run = config["runs"][wildcards.run]
    return expand("mapped-{mapper}/{dataset}.{tissue}.{ref}.bam",
                  dataset=run["dataset"],
                  tissue=tissues,
                  ref=run["ref"],
                  mapper=run["mapper"])


def get_ref(wildcards):
    ref = config["runs"][wildcards.run]["ref"]
    return config["ref"][ref]["fasta"]


wildcard_constraints:
    caller="|".join(config["caller"])


rule all:
    input:
        expand("adhoc-calls/{caller}/{run}.bcf", caller=config["caller"],
                                                 run=config["runs"])


rule test:
    input:
        expand("adhoc-calls/{caller}/{run}.bcf", caller=config["caller"],
                                                 run="test")


include: "rules/mapping.smk"
include: "rules/delly.smk"
include: "rules/pindel.smk"
include: "rules/adhoc.smk"
