from svgutils.compose import *

plots = snakemake.input

Figure("22cm", "6cm",
       SVG(plots[0]),
       SVG(plots[1]).move(100, 0),
       SVG(plots[2]).move(200, 0)
).save(snakemake.output[0])
