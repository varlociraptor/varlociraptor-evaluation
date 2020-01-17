rule get_region_bed:
    input:
        ref=get_ref
    output:
        "neusomatic/{run}.region.bed"
    conda:
        "../envs/pyfaidx.yaml"
    shell:
        "faidx --transform bed {input} > {output}"


rule neusomatic:
    input:
        ref=get_ref,
        fai=lambda w: get_ref(w) + ".fai",
        bed="neusomatic/{run}.region.bed",
        bams=get_bams,
        bais=get_bais,
    output:
        workdir=directory("neusomatic/{run}"),
        vcf="neusomatic/{run}.all.vcf"
    log:
        "logs/neusomatic/{run}.log"
    benchmark:
        "benchmarks/neusomatic/{run}.tsv"
    singularity:
        "docker://msahraeian/neusomatic:0.2.1"
    threads: 20
    shell:
        """
        (preprocess.py --mode call --reference {input.ref} \
                      --region_bed {input.bed} --tumor_bam {input.bams[0]} \
                      --normal_bam {input.bams[1]} --work {output.workdir} \
                      --min_mapq 10 --num_threads {threads} \
                      --scan_alignments_binary /opt/neusomatic/neusomatic/bin/scan_alignments

        call.py --candidates_tsv {output.workdir}/dataset/*/candidates*.tsv \
                --reference {input.ref} --out {output.workdir} \
                --checkpoint /opt/neusomatic/neusomatic/models/NeuSomatic_v0.1.4_standalone_SEQC-WGS-Spike.pth \
                --num_threads {threads} \
                --batch_size 100

        python `which postprocess.py` --reference {input.ref} --tumor_bam {input.bams[0]} \
                       --pred_vcf {output.workdir}/pred.vcf \
                       --candidates_vcf {output.workdir}/work_tumor/filtered_candidates.vcf \
                       --output_vcf {output.vcf} \
                       --work {output.workdir}) 2> {log}
        """


rule fix_neusomatic:
    input:
        header=lambda w: "resources/neusomatic.{}.header.txt".format(config["runs"][w.run]["ref"]),
        vcf="neusomatic/{run}.all.vcf"
    output:
        "default-neusomatic/{run}.all.bcf"
    conda:
        "../envs/tools.yaml"
    shell:
        "bcftools annotate -Ob -o {output} -h {input.header} {input.vcf}"
 


ruleorder: neusomatic_adhoc > adhoc_filter


rule neusomatic_adhoc:
    input:
        "default-neusomatic/{run}.all.bcf"
    output:
        "adhoc-neusomatic/{run}.all.bcf"
    params:
        "-f PASS -Ob"
    wrapper:
        "0.19.3/bio/bcftools/view"
