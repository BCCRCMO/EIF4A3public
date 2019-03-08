#!/usr/bin/Rscript

library(WGCNA)
library(Biobase)
library(ggplot2)
library(tools)


args <- commandArgs(trailingOnly=TRUE)

d <- read.table(args[1], stringsAsFactors=FALSE, header=TRUE)
df <- d[, c(2:ncol(d))]
numSamples <- ncol(df)
rownames(df) <- make.names(d[, 1], unique=TRUE)
df[, "max"] <- apply(d[, 2:(numSamples+1)], 1, max)
df[, "min"] <- apply(d[, 2:(numSamples+1)], 1, min)
df <- df[, 1:numSamples]
print(dim(df))

datExpr <- as.data.frame(t(df))
gsg <- goodSamplesGenes(datExpr, verbose=3)
print(gsg$allOK)

if (!gsg$allOK)
{
  if(sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse=", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse=", ")))
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}

net = blockwiseModules(datExpr, power=as.numeric(args[2]), networkType="signed", TOMType="signed", maxPOutliers=1, numericLabels=TRUE, pamRespectsDendro=FALSE, verbose=3)

library(reshape)

b <- t(datExpr)
b <- t(scale(t(b)))
b <- as.data.frame(b)
b$color <- net$colors
b$names <- rownames(b)
print(table(b$color))

pref <- file_path_sans_ext(args[1])
b_sel <- b[which(b$color != 0), ]
df <- melt.data.frame(b_sel, measure.vars=c("FPKM", paste("FPKM", 1:(numSamples-1), sep='.')))

write.table(unique(df[, c(2, 1)]), quote=FALSE, row.names=FALSE, file=paste(pref, ".txt", sep=''))

counts <- c()
for(i in 1:nrow(df)){
  counts <- c(counts, length(which(df$color == df[i, "color"]))/numSamples)
}

t <- length(unique(df$color))
cs <- rep(0, t)
for(i in 1:length(df$color)){
  cs[df$color[i]] <- counts[i]
}

l <- paste("Group=", as.character(1:t), ", N=", as.character(cs), sep='')
df$color_2 <- factor(paste("Group=", as.character(df$color), ", N=", as.character(counts), sep=''), levels=l)

t <- 0
while (length(which(as.numeric(df$color) <= t)) < nrow(df)*0.8){
      t <- t+1
}
d <- df[which(as.numeric(df$color) <= t), ]
saveRDS(df, file=paste(args[1], ".rds", sep=''))


num_cols <- ceiling(sqrt(t))
num_rows <- ceiling(t/num_cols)

pdf(paste(args[1], ".pdf", sep=''), width=4*num_cols, height = 4*num_rows)
ggplot(d, aes(variable, value, group=names)) + geom_line(size=0.01, alpha=0.1) + stat_summary(fun.y=mean, colour="#08519c", geom="line", size = 3, group=1) + facet_wrap(~color_2, nrow=num_rows, ncol=num_cols) + theme_bw() + ylab(expression(paste(psi, " (normalized)"))) + xlab(expression(paste("Dose (", mu, "u)"))) + scale_x_discrete(labels = c(0, rep(c(0.5, 2, 5, 10, 20), 1))) + theme(text = element_text(size = 20), axis.text.x=element_text(angle=90, vjust=0.4), panel.border = element_rect(colour="black", size = 1))
dev.off()
