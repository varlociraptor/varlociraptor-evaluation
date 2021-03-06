library(UpSetR)
library(plyr)
library(ggplot2)
library(grid)
library(rlist)

calls <- read.table(snakemake@input[[1]], row.names = NULL, sep = "\t", header = TRUE, check.names = FALSE)
#print(dim(calls))
calls <- calls[calls$max_case_af >= 0.1, ]
#calls <- calls[calls$max_prob_somatic_tumor < 0.07]_
#calls$max_case_af_log10 <- log10(calls$max_case_af)

queries <- list()
for(ds in snakemake@params[["datasets"]]) {
    queries <- list.append(queries, list(query=intersects, params=list(ds)))
}


svg(file = snakemake@output[[1]])
if(!empty(calls)) {
    upset(calls, 
          order.by = "freq", 
          sets = snakemake@params[["datasets"]],
          queries = queries
    )
}
dev.off()
