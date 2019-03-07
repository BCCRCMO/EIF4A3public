#!/share/lustre/backup/amazloomian/Projects/install/R/R-3.1.2/bin/Rscript

library(edgeR)

args <- commandArgs(trailingOnly = TRUE)

expressions <- read.table(args[1], stringsAsFactors=FALSE, header=TRUE)
d <- expressions[, 3:ncol(expressions)]
normalized = t(t(d)*calcNormFactors(d, method="upperquartile"))
d <- cbind(expressions[, 1:2], normalized)

write.table(d, quote=FALSE, row.names=FALSE)
