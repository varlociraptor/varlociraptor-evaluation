from pybedtools import BedTool


assert len(snakemake.input.bams) == 4

MINCOV = 10

exons = BedTool(snakemake.input.exons)

coverage = exons.multi_bam_coverage(bams=snakemake.input.bams)


def filter(exon):
    return all(cov / len(exon) >= MINCOV for cov in exon[-4:])
        

coverage.filter(filter).saveas(snakemake.output.bed)
