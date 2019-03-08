#!/usr/bin/Rscript

library(ggplot2)

d <- read.table("cluster_events")
colnames(d) <- c("Library", "Type", "Direction")

pdf("event_counts_in_clusters.pdf")
ggplot(d, aes(Direction)) + geom_bar(aes(fill=Type), position="fill", color="black") + facet_wrap(~Library) + xlab("\nDirection") + ylab("Ratio\n") + theme_bw() + theme(text = element_text(size = 16), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) + scale_fill_brewer(palette = "Set3")
dev.off()
