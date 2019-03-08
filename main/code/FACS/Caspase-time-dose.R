
### R script used
wd="/Users/dyap/Documents/Collboration_Projects/EIF4A3/FACS"
setwd(wd)
T595<- read.csv("T-595-48.txt", header=FALSE)
T598<- read.csv("T-598-48.txt", header=FALSE)

# Take the first row as the header, to keep it in numerical format
colnames(T595)<-round(T595[1,], digits = 2)
T595 = T595[-1,]
colnames(T598)<-round(T598[1,], digits = 2)
T598 = T598[-1,]

# Order the means by colnames
T595_means <- (colMeans(T595, na.rm = TRUE))
T595.order <- order (T595_means)
ordered.names <- names(T595_means)[T595.order]

T598.order <- match(ordered.names, names(T598))
T598_means <- (colMeans(T598, na.rm = TRUE))

pdf(file = "48hr_Caspase.pdf", width =10, height = 7)

par(mar = c(5,5,5,5))

boxplot(T595 [, order (T595_means)], boxwex = 0.25, at= 1:10 - 0.2, horizontal = FALSE, outline = FALSE, cex=3, cex.axis=1.2
, las=1, xlim = c(0, 11), ylim = c(0, 35000))
points(T595_means [ order (T595_means)], pch = 18, col = "blue", lwd = 2, cex=1.3)
mtext("Caspase 3/7 Activity", line=4, cex=1.3, side=2, col="black")
text(5.5,34000, "48 hr", cex=2.0, col="black")

stripchart(T595 [, order (T595_means)], vertical = TRUE, 
    method = "jitter", at= 1:10 - 0.2, add = TRUE, pch = 18, col = 'blue')

boxplot(T598 [, order (T598.order)], boxwex = 0.25, at = 1:10 + 0.2, horizontal = FALSE, outline = FALSE, xaxt='n', ann=FALSE,
yaxt='n', xlab="", ylab="", add = TRUE, xlim = c(0, 11), ylim = c(0, 35000))
points(T598_means [ order (T598_means)], pch = 20, col = "red", lwd = 2)

stripchart(T598 [, order (T598.order)], vertical = TRUE,
    method = "jitter", at = 1:10 + 0.2, add = TRUE, pch = 20, col = 'red', cex=1.3)
    
# fit a loess line
#loess_fit <- loess(T595_means ~ T595.order, T595)
#lines(T595.order, predict(loess_fit), col = "blue")

#loess_fit <- loess(T598_means ~ T598.order, T598)
#lines(T598.order, predict(loess_fit), col = "red")
    
# Legend    
    legend(0.5, 34500, c("T-595", "T-598"),
       col = c("blue", "red"), pch = c(18,20), cex=1.3)


dev.off()
    
    
#############################################