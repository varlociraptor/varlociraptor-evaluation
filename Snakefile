import pandas as pd

from snakemake.remote import FTP

ftp = FTP.RemoteProvider()


configfile: "config.yaml"


gdc_manifest = pd.read_table(config["gdc-manifest"])
units = pd.DataFrame({
    "sample": gdc_manifest.filename.str.slice(5, 33),
    "path": expand("gdc-data/{bam}", bam=gdc_manifest.filename)})
units.index = gdc_manifest.id


tissues = ["tumor", "normal"]
# pattern for accessing bam files given a run
bams = expand("mapped/{{run}}.{tissue}.bam", tissue=tissues)


rule all:
    input:
        expand("mapped/{run}.", unit=units.index)


include: "rules/mapping.smk"
include: "rules/delly.smk"
include: "rules/pindel.smk"
