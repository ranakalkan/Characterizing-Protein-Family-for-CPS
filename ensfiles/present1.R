setwd("/home/fatih/Documents/ENS210/project/fastas")

library('seqinr')
library('alignfigR')
library('msa')
'''
###SEQINR###
seqData <- read.alignment('more_simplified.fas', format='fasta')
dist <- seqinr::dist.alignment(seqData)
scatter.smooth(dist)

matrix <- as.matrix(dist.alignment(seqData, matrix = "identity" ))
heatmap(matrix,margins=c(16, 8))
###SEQINR###

###ALIGNFIGR###
library(grid)
library(ggplot2)
figrdata = read.fasta('more_simplified.fas')
alignfigR::plot_alignment(figrdata)
grob <- grobTree(textGrob("CSP1", x=0.9,  y=0.5, hjust=0,
                          gp=gpar(col="black", fontsize=8, fontface="italic")))
plot_alignment(figrdata) +  annotation_custom(grob) 


plot_alignment(figrdata) + annotate("text", x = -120, y = 11, label="CPS2") + annotate("text", x = -120, y = 15, label="CPS3") + annotate("text", x = -120, y = 26, label="CPS1") + annotate("segment", x = 0, ,y= 11, ,xend = 2355, yend = 11, colour = "black") + annotate("segment", x = 0, ,y= 15, ,xend = 2355, yend = 15, colour = "black") + annotate("segment", x = 0, ,y= 26, ,xend = 2355, yend = 26, colour = "black")
'''
###ALIGNFIGR###

###MSA###
library(msa)
system.file("tex", "texshade.sty", package="msa")

CPSI <- readAAStringSet("fasta_110_CPSI.fas")
CPSII <- readAAStringSet("fasta_110_CPSII.fas")
CPSIII <- readAAStringSet("fasta_4_CPSIII.fas")
CPS_I_II <- readAAStringSet("fasta_110_CPS_I_II.fas")
CPS_I_III <- readAAStringSet("fasta_110_CPS_I_III.fas")
CPS_II_III <- readAAStringSet("fasta_110_CPS_II_III.fas")
CPS_I_II_III <- readAAStringSet("fasta_110_CPS_I_II_III.fas")

I_aligned <- msa(CPSI)
II_aligned <- msa(CPSII)
III_aligned <- msa(CPSIII)

I_II_aligned <- msa(CPS_I_II)
I_III_aligned <- msa(CPS_I_III)
II_III_aligned <- msa(CPS_II_III)
I_II_III_aligned <- msa(CPS_I_II_III)

I_consensus <- msa::msaConsensusSequence(I_aligned)
II_consensus <- msa::msaConsensusSequence(II_aligned)
I_II_consensus <- msa::msaConsensusSequence(I_II_aligned)


msaData

msaPrettyPrint(msaData, y=c(100, 213), output="pdf",
               showNames="none", showLogo="top", askForOverwrite=FALSE)
  
install.packages("ape")
library(ape)
#d <- dist.alignment(hemoAln2, "identity")
hemoTree <- nj(dist)
plot(hemoTree, main="Phylogenetic Tree of CSP Sequences")
data(BLOSUM62)

conservation_CPS_I = msa::msaConservationScore(I_aligned, BLOSUM62)
conservation_CPS_II = msa::msaConservationScore(II_aligned, BLOSUM62)
conservation_CPS_III = msa::msaConservationScore(III_aligned, BLOSUM62)
conservation_CPS_I_II = msa::msaConservationScore(I_II_aligned, BLOSUM62)
conservation_CPS_I_III = msa::msaConservationScore(I_III_aligned, BLOSUM62)
conservation_CPS_II_III = msa::msaConservationScore(II_III_aligned, BLOSUM62)
conservation_CPS_I_II_III = msa::msaConservationScore(I_II_III_aligned, BLOSUM62)

#par(mfrow=c(2,3)) # all plots on one page

attach(mtcars)
layout(matrix(c(1,2,3,4,5,6,7,7,7), 3, 3, byrow = TRUE))

plot(conservation_CPS_I, type="n")
lines(conservation_CPS_I, type="h")

plot(conservation_CPS_II, type="n")
lines(conservation_CPS_II, type="h")

plot(conservation_CPS_III, type="n")
lines(conservation_CPS_III, type="h")

plot(conservation_CPS_I_II, type="n")
lines(conservation_CPS_I_II, type="h")

plot(conservation_CPS_I_III, type="n")
lines(conservation_CPS_I_III, type="h")

plot(conservation_CPS_II_III, type="n")
lines(conservation_CPS_II_III, type="h")

plot(conservation_CPS_I_II_III, type="n")
lines(conservation_CPS_I_II_III, type="h")

##1 Row'da consensus'ların msa'ları

conses

###MSA###

