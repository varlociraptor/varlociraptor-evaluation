from cyvcf2 import VCF, Writer
import numpy as np


def subclone_vaf(gt):
    """Calculate subclone allele frequency"""
    if np.all(gt[:2] == [1, 1]):
        return 1.0
    elif (np.all(gt[:2] == [0, 1]) or np.all(gt[:2] == [1, 0]) or
          np.all(gt[:2] == [-1, 1]) or np.all(gt[:2] == [1, -1])):
        return 0.5
    else:
        return 0.0


# Reader
vcf_in = VCF(snakemake.input[0])

# Setup subclone information
subclones = ["Som{}".format(i) for i in range(1, 5)]
fractions = [1/3, 1/3, 1/4, 1/12]


# Prepare writer
vcf_in.add_info_to_header({"ID": "TAF",
                           "Number": "1",
                           "Description": "True tumor allele frequency",
                           "Type": "Float"})
vcf_in.add_info_to_header({"ID": "NAF",
                           "Number": "1",
                           "Description": "True normal allele frequency",
                           "Type": "Float"})
bcf_out = Writer(snakemake.output[0], vcf_in)

for rec in vcf_in:
    if len(rec.ALT) > 1:
        raise ValueError("multiallelic sites are not supported at the moment")

    try:
        # get VAFs from VCF
        tumor_vaf = rec.INFO["TAF"]
        normal_vaf = rec.INFO["NAF"]
    except KeyError:
        # calculate VAFs
        subclone_idx = [vcf_in.samples.index(s) for s in subclones]
        control_idx = vcf_in.samples.index("Control")

        tumor_vaf = sum(fraction * subclone_vaf(rec.genotypes[idx])
                        for idx, fraction in zip(subclone_idx, fractions))
        normal_vaf = subclone_vaf(rec.genotypes[control_idx])

        rec.INFO["TAF"] = tumor_vaf
        rec.INFO["NAF"] = normal_vaf
        
    # only keep somatic variants
    if normal_vaf == 0.0 and tumor_vaf > 0.0:
        bcf_out.write_record(rec)

bcf_out.close()
