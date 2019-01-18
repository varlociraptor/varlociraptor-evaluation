library(UpSetR)
library(plyr)
library(ggplot2)
library(grid)
library(rlist)

calls <- read.table(snakemake@input[[1]], row.names = NULL, sep = "\t", header = TRUE, check.names = FALSE)
#calls <- calls[calls$max_case_af >= 0.2, ]
calls$max_case_af_log10 <- log10(calls$max_case_af)

queries <- list()
for(ds in snakemake@params[["datasets"]]) {
    queries <- list.append(queries, list(query=intersects, params=list(ds)))
}
print(queries)


svg(file = snakemake@output[[1]])
upset(calls, 
      order.by = "freq", 
      sets = snakemake@params[["datasets"]],
      attribute.plots = list(
          gridrows = 20, ncols = 2,
          plots = list(
              list(plot=histogram, x="max_case_af_log10", queries=TRUE),
              list(plot=scatter_plot, x="max_case_af_log10", y="max_prob_somatic_tumor", queries=TRUE)
          )
      ),
      queries = queries
)
dev.off()
