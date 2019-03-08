#!/usr/bin/Rscript

library(ggplot2)
library(grid)

give.n <- function(x){
  return(data.frame(y = as.numeric(quantile(x, 0.95))+0.1, label = paste("N=", length(x), sep='')))
}


d <- read.table("all_events", stringsAsFactors=FALSE)

pdf("Dose_histogram_all.pdf", width=12, height=6)
ggplot(d, aes(factor(V6))) + geom_bar() + facet_wrap(~V1, nrow=2) + theme_bw() + theme(text = element_text(size = 16), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), panel.border = element_rect(colour="black", size = 1)) + scale_x_discrete(labels = c("0.5", "2", "5", "10", "20")) + xlab(expression(paste("Dose (", mu, "u)"))) + ylab("Count\n")
dev.off()

df <- within(d, { count <- ave(V1, V1, V2, V6, FUN=function(x) length(x))})[, c(6, 2, 8, 1)]
df <- unique(df)
colnames(df) <- c("V1", "V2", "V3", "V4")
df$V3 <- as.numeric(df$V3)

pdf("Dose_type_histogram_all.pdf", width=12, height=6)
ggplot(df, aes(V1, V3, fill=V2)) + geom_line(aes(color=factor(V2)), size=1.5, linetype = 1) + geom_point(size=3, shape=24) + theme_bw() + theme(panel.spacing = unit(1, "lines"), text = element_text(size = 16), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), panel.border = element_rect(colour="black", size = 1)) + ylab("Count\n") + scale_x_continuous(name = expression(paste("Drug level (", mu, "u)")), breaks = c(0, 1, 2, 3, 4, 5), labels = c("0", "0.5", "2", "5", "10", "20")) + scale_fill_discrete(name = "Event type\n") + scale_color_discrete(guide=FALSE) + facet_wrap(~V4, nrow=2)
dev.off()

exps <- unique(d[, 1])

for(i in 1:length(exps)){
  df <- d[which(d[, 1] == exps[i]), ]
  pdf(paste(exps[i], "_diff_boxplot.pdf", sep=''))
  p1 <- ggplot(df, aes(factor(V6), V4)) + geom_boxplot(outlier.shape=NA) + stat_summary(fun.data=give.n, geom = "text", fun.y = median, size = 3, angle=30) + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + facet_wrap(~V2) + scale_x_discrete(labels = c("0.5", "2", "5", "10", "20")) + xlab(expression(paste("Dose (", mu, "u)"))) + ylab(expression(paste(Delta, " (", psi, ")"))) + ggtitle(exps[i])
  print(p1)
  dev.off()
}

d <- read.table("diff_events_per_drug", stringsAsFactors=FALSE, header=FALSE)
exps <- unique(d[, 1])

for(i in 1:length(exps)){
  df <- d[which(d[, 1] == exps[i]), ]
  pdf(paste(exps[i], "_response_boxplot.pdf", sep=''))
  p1 <- ggplot(df, aes(factor(V3), V4)) + geom_boxplot(outlier.shape=NA) + stat_summary(fun.data=give.n, geom = "text", fun.y = median, size = 3, angle=30) + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + facet_wrap(~V2) + scale_x_discrete(labels = c("0.5", "2", "5", "10", "20")) + xlab(expression(paste("Dose (", mu, "u)"))) + ggtitle(paste(exps[i], sep='')) + ylab(expression(paste("|", psi[i], "|  -  |", psi[i-1], "|")))
  print(p1)
  dev.off()
}
