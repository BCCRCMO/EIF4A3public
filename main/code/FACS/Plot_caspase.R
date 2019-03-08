### Stress granules graphs
sg24df <- read.table(file = "../../data/Stress_Granules/Hela_analysis_24hr.csv", col.names = c("Grp", "SGprp", "Trt"), sep = ",", stringsAsFactors = FALSE)

sg24df$Trtf <- factor(sg24df$Trt, levels = c("UT", "T-595", "T-598"))
sg24df$SGpct <- sg24df$SGprp*100

pdf(file = "SG_HeLa_24hr.pdf",  width = 4, height = 6, useDingbats = FALSE)
plot(sg24df$Trtf, sg24df$SGpct,
     ylab = "Cells displaying SGs(%)",
     xlab = "24hr inhibitor (last 1hr with Ars)",
     col = c("grey40", "white", "red"),
     ylim = c(0, 119), yaxt = "n", lwd = 1.3
     )
axis(side=2, at = c(20, 40, 60, 80, 100))
points(sg24df$Trtf, sg24df$SGpct, pch=19)
segments(1, 102, 1, 105, lwd = 3)
segments(1, 105, 2, 105, lwd = 3)
segments(2, 105, 2, 65, lwd = 3)
text(x=1.5, y = 108, labels="p = 0.023")
segments(1, 112, 1, 115, lwd = 3)
segments(1, 115, 3, 115, lwd = 3)
segments(3, 115, 3, 100, lwd = 3)
text(x=2.5, y = 118, labels="p = 0.27")
dev.off()

###

sg04df <- read.table(file = "../../data/Stress_Granules/Hela_analysis_04hr.csv", col.names = c("Grp", "SGavg", "Trt"), sep = ",", stringsAsFactors = FALSE)

sg04df$Trtf <- factor(sg04df$Trt, levels = c("UT", "T-595", "T-598"))


pdf(file = "SG_HeLa_04hr.pdf",  width = 4, height = 6, useDingbats = FALSE)
plot(sg04df$Trtf, sg04df$SGavg,
     ylab = "Cells displaying SGs(%)",
     xlab = "4hr inhibitor (last 1hr with Ars)",
     col = c("grey40", "white", "red"),
     ylim = c(0, 13), yaxt = "n", lwd = 1.3
     )
axis(side=2, at = c(1:12))
points(sg04df$Trtf, sg04df$SGavg, pch=19)
segments(1, 10.2, 1, 10.5, lwd = 3)
segments(1, 10.5, 2, 10.5, lwd = 3)
segments(2, 10.5, 2, 6.5, lwd = 3)
text(x=1.5, y = 10.8, labels="p = 0.044")
segments(1, 11.2, 1, 12.5, lwd = 3)
segments(1, 12.5, 3, 12.5, lwd = 3)
segments(3, 12.5, 3, 11.5, lwd = 3)
text(x=2.5, y = 12.8, labels="p = 0.65")
dev.off()
