#!/usr/bin/Rscript

library(googleVis)
library(plyr)

f <- function(From, To){
  a <- scan(as.character(From), character())
  b <- scan(as.character(To), character())
  return(length(intersect(a, b)))
}

##################################################
# Cell line vs. Cell line
##################################################

Hela <- c("Hela_T_202.up", "Hela_T_202.down", "Hela_T_595.up", "Hela_T_595.down")
HCT116 <- c("HCT116_T_202.up", "HCT116_T_202.down", "HCT116_T_595.up", "HCT116_T_595.down")
cell_line_overlaps = expand.grid(From = Hela, To = HCT116)
cell_line_overlaps = mdply(cell_line_overlaps, f)
colnames(cell_line_overlaps) <- c("From", "To", "Weight")

Sankey <- gvisSankey(cell_line_overlaps, from="From", to="To", weight="Weight",
       	  				 list(
                       sankey="{
                           node: { colors: ['#d73027', '#d73027', '#4575b4', '#d73027', '#4575b4', '#4575b4', '#d73027', '#4575b4'], width: 30,
                            label: { color: '#000', fontSize: 12, fontName: 'Times-Roman'} }}"))


print(Sankey, file="cell_lines_sankey.html")
system("wkhtmltopdf --enable-plugins --javascript-delay 10000  cell_lines_sankey.html cell_lines_sankey.pdf")

##################################################
# Compound vs. compound
##################################################

T_595 <- c("Hela_T_595.up", "Hela_T_595.down", "HCT116_T_595.up", "HCT116_T_595.down")
T_202 <- c("Hela_T_202.up", "Hela_T_202.down", "HCT116_T_202.up", "HCT116_T_202.down")
compound_overlaps = expand.grid(From = T_202, To = T_595)
compound_overlaps = mdply(compound_overlaps, f)
colnames(compound_overlaps) <- c("From", "To", "Weight")

Sankey <- gvisSankey(compound_overlaps, from="From", to="To", weight="Weight",
       	  				 list(
                       sankey="{
                           node: { colors: ['#d73027', '#d73027', '#4575b4', '#d73027', '#4575b4', '#4575b4', '#d73027', '#4575b4'], width: 30,
                            label: { color: '#000', fontSize: 12, fontName: 'Times-Roman'} }}"))


print(Sankey, file="compounds_sankey.html")
system("wkhtmltopdf --enable-plugins --javascript-delay 10000  compounds_sankey.html compound_sankey.pdf")
