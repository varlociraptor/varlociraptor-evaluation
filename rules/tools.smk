rule svg2pdf:
    input:
        "{prefix}.svg"
    output:
        "{prefix}.pdf"
    conda:
        "../envs/cairosvg.yaml"
    shell:
        "cairosvg {input} -o {output}"
