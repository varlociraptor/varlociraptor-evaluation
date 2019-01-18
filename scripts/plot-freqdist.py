from pysam import VariantFile
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import islice
import pandas as pd


sns.set_style("ticks")


calls = VariantFile(snakemake.input[0])


def freq(sample):
    ref, alt = sample.get("AD")
    if alt == 0:
        return 0
    return alt / (ref + alt)


normal_freqs = []
tumor_freqs = []
vartypes = []

print("collecting calls")
for record in calls:
    if not record.info.get("SHARED"):
        continue
    normal_freqs.append(freq(record.samples.get("normal")))
    tumor_freqs.append(freq(record.samples.get("tumor")))
    vartypes.append(record.info.get("TYPE"))

d = pd.DataFrame({"normal": normal_freqs, "tumor": tumor_freqs, "type": vartypes})

print("subsampling")
# sample d to not get overwhelmed
d = d.sample(10000, random_state=245746)

print("plotting and density estimation")
g = sns.FacetGrid(col="type", data=d)
g.map(sns.kdeplot, "normal", "tumor", shade=True, clip=(0, 1), shade_lowest=False)
g.map(sns.scatterplot, "normal", "tumor", size=1, alpha=0.5, marker=".")

plt.savefig(snakemake.output[0], bbox_inches="tight")
