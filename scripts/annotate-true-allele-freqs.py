from cyvcf2 import VCF, Writer


def subclone_vaf(gt):
    """Calculate subclone allele frequency"""
    if np.all(gt[:2] == [1, 1]):
        return 1.0
    elif np.all(gt[:2] == [0, 1]) or np.all(gt[:2] == [1, 0]):
        return 0.5
    else:
        return 0.0


# Reader
vcf_in = VCF(snakemake.input.vcf)

# Setup subclone information
subclones = ["Som{}".format(i) for i in range(1, 5)]
fractions = [1/3, 1/3, 1/4, 1/12]
subclone_idx = [vcf_in.samples.index(s) for s in subclones]


# Prepare writer
bcf_out = Writer(snakemake.output[0], bcf_in)
bcf_out.add_info_to_header({"ID": "AF",
                            "Number": "1",
                            "Description": "True tumor allele frequency",
                            "Type": "Float"})

for rec in vcf_in:
    if len(rec.ALT) > 1:
        raise ValueError("multiallelic sites are not supported at the moment")

    # calculate AF
    vaf = sum(fraction * subclone_vaf(rec.genotypes[idx])
              for idx, fraction in zip(subclone_idx, fractions))

    rec.INFO["AF"] = vaf
    bcf_out.write_record(rec)

bcf_out.close()