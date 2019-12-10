setwd("/home/fatih/Desktop/ensfiles")

library('seqinr')
library('alignfigR')
library('msa')


###SEQINR###
seqData <- read.alignment('simplified.fas', format='fasta')
dist <- seqinr::dist.alignment(seqData)
scatter.smooth(dist)

matrix <- as.matrix(dist.alignment(seqData, matrix = "identity" ))
heatmap(matrix,margins=c(16, 8))
###SEQINR###

###ALIGNFIGR###
library(grid)
library(ggplot2)
figrdata = read.fasta('simplified.fas')
alignfigR::plot_alignment(figrdata)
grob <- grobTree(textGrob("CSP1", x=0.9,  y=0.5, hjust=0,
                          gp=gpar(col="black", fontsize=8, fontface="italic")))
plot_alignment(figrdata) +  annotation_custom(grob) 


plot_alignment(figrdata) + annotate("text", x = -120, y = 11, label="CPS2") + annotate("text", x = -120, y = 15, label="CPS3") + annotate("text", x = -120, y = 26, label="CPS1") + annotate("segment", x = 0, ,y= 11, ,xend = 2355, yend = 11, colour = "black") + annotate("segment", x = 0, ,y= 15, ,xend = 2355, yend = 15, colour = "black") + annotate("segment", x = 0, ,y= 26, ,xend = 2355, yend = 26, colour = "black")
###ALIGNFIGR###

###MSA###
library(msa)
system.file("tex", "texshade.sty", package="msa")

mySequences <- readAAStringSet("simplified.fas")
#mySequences <- readAAStringSet("AllCPSsUnaligned.fas")
msaData <- msa(mySequences)

msaPrettyPrint(msaData, y=c(100, 213), output="pdf",
               showNames="none", showLogo="top", askForOverwrite=FALSE)

install.packages("ape")
library(ape)
#d <- dist.alignment(hemoAln2, "identity")
hemoTree <- nj(dist)
plot(hemoTree, main="Phylogenetic Tree of CSP Sequences")
###MSA###


