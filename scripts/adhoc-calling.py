from cyvcf2 import VCF, Writer

def get_sample_name(tissue):
    ds = snakemake.config["runs"][snakemake.wildcards.run]["dataset"]
    return snakemake.config["datasets"][ds][tissue]["name"]


tumor, normal = map(get_sample_name, ["tumor", "normal"])


bcf_in = VCF(snakemake.input.vcf)
bcf_out = Writer(snakemake.output[0], bcf_in)


for rec in bcf_in:
    if rec.FILTER:
        continue
    gt = rec.genotypes
    tumor_gt = gt[0, :2]
    normal_gt = gt[1, :2]
    if (np.any(tumor_gt) and
        not np.any(normal_gt) and
        not np.any(np.isnan(normal_gt))):
        # somatic variant
        bcf_out.write_record(rec)
