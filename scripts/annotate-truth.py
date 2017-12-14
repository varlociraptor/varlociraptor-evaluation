from cyvcf2 import VCF, Writer
import numpy as np


def subclone_vaf(gt):
    """Calculate subclone allele frequency"""
    if np.all(gt[:2] == [1, 1]):
        return 1.0
    elif np.all(gt[:2] == [0, 1]) or np.all(gt[:2] == [1, 0]):
        return 0.5
    else:
        return 0.0


# Reader
vcf_in = VCF(snakemake.input[0])

# Setup subclone information
subclones = ["Som{}".format(i) for i in range(1, 5)]
fractions = [1/3, 1/3, 1/4, 1/12]
subclone_idx = [vcf_in.samples.index(s) for s in subclones]
control_idx = vcf_in.samples.index("Control")


# Prepare writer
vcf_in.add_info_to_header({"ID": "AF",
                            "Number": "1",
                            "Description": "True tumor allele frequency",
                            "Type": "Float"})
vcf_in.add_info_to_header({"ID": "SOMATIC",
                            "Number": "0",
                            "Description": "Somatic mutation",
                            "Type": "Flag"})
bcf_out = Writer(snakemake.output[0], vcf_in)

for rec in vcf_in:
    if len(rec.ALT) > 1:
        raise ValueError("multiallelic sites are not supported at the moment")

    # calculate AF
    vaf = sum(fraction * subclone_vaf(rec.genotypes[idx])
              for idx, fraction in zip(subclone_idx, fractions))

    print(vaf)
    rec.INFO["AF"] = vaf
    rec.INFO["SOMATIC"] = (subclone_vaf(rec.genotypes[control_idx]) == 0.0 and
                           vaf > 0.0)
    bcf_out.write_record(rec)

bcf_out.close()
