#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
d <- read.table(args[1], stringsAsFactors = FALSE)
# print(head(d))
d["p_value"] <- apply(d, 1, function(x) phyper(max(as.numeric(x[2]), 0), as.numeric(x[4]), as.numeric(x[5])-as.numeric(x[4]), as.numeric(x[3]), lower.tail=FALSE))
d["FDR"] <- p.adjust(d[["p_value"]], "BH")
d <- d[order(d$FDR), ]
cat("ID p_value FDR\n")
cat(paste(d[["V1"]], d[["p_value"]], d[["FDR"]], sep=' '), sep='\n')
