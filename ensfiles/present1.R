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
library(ggplot2)
system.file("tex", "texshade.sty", package="msa")

CPSI <- readAAStringSet("fasta_110_CPSI.fas")
CPSII <- readAAStringSet("fasta_110_CPSII.fas")
CPSIII <- readAAStringSet("fasta_4_CPSIII.fas")
CPS_I_II <- readAAStringSet("fasta_110_CPS_I_II.fas")
CPS_I_III <- readAAStringSet("fasta_110_CPS_I_III.fas")
CPS_II_III <- readAAStringSet("fasta_110_CPS_II_III.fas")
CPS_I_II_III <- readAAStringSet("fasta_110_CPS_I_II_III.fas")

consensus <- readAAStringSet("conses.fas")

I_aligned <- msa(CPSI)
II_aligned <- msa(CPSII)
III_aligned <- msa(CPSIII)

I_II_aligned <- msa(CPS_I_II)
I_III_aligned <- msa(CPS_I_III)
II_III_aligned <- msa(CPS_II_III)
I_II_III_aligned <- msa(CPS_I_II_III)

cons_aligned <- msa(consensus)

writeXStringSet(unmasked(I_aligned), filepath = "I_aligned.fas")
writeXStringSet(unmasked(II_aligned), filepath = "II_aligned.fas")
writeXStringSet(unmasked(III_aligned), filepath = "III_aligned.fas")

writeXStringSet(unmasked(I_II_aligned), filepath = "I_II_aligned.fas")
writeXStringSet(unmasked(I_III_aligned), filepath = "I_III_aligned.fas")
writeXStringSet(unmasked(II_III_aligned), filepath = "II_III_aligned.fas")

writeXStringSet(unmasked(I_II_III_aligned), filepath = "I_II_III_aligned.fas")

cat(">CPSI_consensus\n", file="conses.fas",append=TRUE)
cat(I_consensus, file="conses.fas",append=TRUE)
cat("\n", file="conses.fas",append=TRUE)

cat(">CPSII_consensus\n", file="conses.fas",append=TRUE)
cat(II_consensus, file="conses.fas",append=TRUE)
cat("\n", file="conses.fas",append=TRUE)

cat(">CPSIII_consensus\n", file="conses.fas",append=TRUE)
cat(III_consensus, file="conses.fas",append=TRUE)
cat("\n", file="conses.fas",append=TRUE)



I_consensus <- msa::msaConsensusSequence(I_aligned)
II_consensus <- msa::msaConsensusSequence(II_aligned)
III_consensus <- msa::msaConsensusSequence(III_aligned)
I_II_consensus <- msa::msaConsensusSequence(I_II_aligned)
I_III_consensus <- msa::msaConsensusSequence(I_III_aligned)
II_III_consensus <- msa::msaConsensusSequence(II_III_aligned)
I_II_III_consensus <- msa::msaConsensusSequence(I_II_III_aligned)


data(BLOSUM62)

conservation_CPS_I = msa::msaConservationScore(I_aligned, BLOSUM62)
conservation_CPS_II = msa::msaConservationScore(II_aligned, BLOSUM62)
conservation_CPS_III = msa::msaConservationScore(III_aligned, BLOSUM62)
conservation_CPS_I_II = msa::msaConservationScore(I_II_aligned, BLOSUM62)
conservation_CPS_I_III = msa::msaConservationScore(I_III_aligned, BLOSUM62)
conservation_CPS_II_III = msa::msaConservationScore(II_III_aligned, BLOSUM62)
conservation_CPS_I_II_III = msa::msaConservationScore(I_II_III_aligned, BLOSUM62)
conservation_consensus = msa::msaConservationScore(cons_aligned, BLOSUM62)

vecI <- vector()
vecI <- c(vecI, 1:length(conservation_CPS_I))

vecII <- vector()
vecII <- c(vecII, 1:length(conservation_CPS_II))

vecIII <- vector()
vecIII <- c(vecIII, 1:length(conservation_CPS_III))

vecI_II <- vector()
vecI_II <- c(vecI_II, 1:length(conservation_CPS_I_II))

vecI_III <- vector()
vecI_III <- c(vecI_III, 1:length(conservation_CPS_I_III))

vecII_III <- vector()
vecII_III <- c(vecII_III, 1:length(conservation_CPS_II_III))

vecI_II_III <- vector()
vecI_II_III <- c(vecI_II_III, 1:length(conservation_CPS_I_II_III))

vec_cons <- vector()
vec_cons <- c(vec_cons, 1:length(conservation_consensus))

install.packages("gridExtra")
library(gridExtra)
library(grid)
library(lattice)
library(cowplot)
install.packages("cowplot")

p1 <- ggplot() + geom_line(aes(vecI,conservation_CPS_I)) + annotate("rect", xmin = 55, xmax = 409, ymin = 0, ymax = max(conservation_CPS_I), alpha = .2, fill = 'blue') + annotate("rect", xmin = 430, xmax = 1502, ymin = 0, ymax = max(conservation_CPS_I), alpha = .2, fill = 'orange') + annotate("rect", xmin = 1371, xmax = 1486, ymin = max(conservation_CPS_I) - max(conservation_CPS_I)*0.9, ymax = max(conservation_CPS_I)*0.9, alpha = .2, fill = 'red')  + annotate("rect", xmin = 850, xmax = 973, ymin = max(conservation_CPS_I) - max(conservation_CPS_I)*0.9, ymax = max(conservation_CPS_I)*0.9, alpha = .2, fill = 'cyan')
p2 <- ggplot() + geom_line(aes(vecII,conservation_CPS_II)) + annotate("rect", xmin = 1, xmax = 355, ymin = 0, ymax = max(conservation_CPS_II), alpha = .2, fill = 'blue') + annotate("rect", xmin = 391, xmax = 1440, ymin = 0, ymax = max(conservation_CPS_II), alpha = .2, fill = 'orange') + annotate("rect", xmin = 1313, xmax = 1438, ymin = max(conservation_CPS_II) - max(conservation_CPS_II)*0.9, ymax = max(conservation_CPS_II)*0.9, alpha = .2, fill = 'red')  + annotate("rect", xmin = 799, xmax = 921, ymin = max(conservation_CPS_II) - max(conservation_CPS_II)*0.9, ymax = max(conservation_CPS_II)*0.9, alpha = .2, fill = 'cyan') + annotate("rect", xmin = 1460, xmax = 1806, ymin = 0, ymax = max(conservation_CPS_II), alpha = .2, fill = 'purple') + annotate("rect", xmin = 1812, xmax = 1907, ymin = 0, ymax = max(conservation_CPS_II), alpha = .2, fill = 'yellow1') + annotate("rect", xmin = 1925, xmax = 2065, ymin = 0, ymax = max(conservation_CPS_II), alpha = .4, fill = 'orchid')# + annotate("rect", xmin = 1920, xmax = 2224, ymin = 0, ymax = max(conservation_CPS_I), alpha = .2, fill = 'green1') 
p3 <- ggplot() + geom_line(aes(vecIII,conservation_CPS_III)) + annotate("rect", xmin = 55, xmax = 409, ymin = 0, ymax = max(conservation_CPS_III), alpha = .2, fill = 'blue') + annotate("rect", xmin = 430, xmax = 1502, ymin = 0, ymax = max(conservation_CPS_III), alpha = .2, fill = 'orange') + annotate("rect", xmin = 1371, xmax = 1486, ymin = max(conservation_CPS_III) - max(conservation_CPS_III)*0.9, ymax = max(conservation_CPS_III)*0.9, alpha = .2, fill = 'red')  + annotate("rect", xmin = 850, xmax = 973, ymin = max(conservation_CPS_III) - max(conservation_CPS_III)*0.9, ymax = max(conservation_CPS_III)*0.9, alpha = .2, fill = 'cyan')
p4 <- ggplot() + geom_line(aes(vecI_II,conservation_CPS_I_II)) + annotate("rect", xmin = 1, xmax = 355, ymin = 0, ymax = max(conservation_CPS_I_II), alpha = .2, fill = 'blue') + annotate("rect", xmin = 391, xmax = 1440, ymin = 0, ymax = max(conservation_CPS_I_II), alpha = .2, fill = 'orange') + annotate("rect", xmin = 1313, xmax = 1438, ymin = max(conservation_CPS_I_II) - max(conservation_CPS_I_II)*0.9, ymax = max(conservation_CPS_I_II)*0.9, alpha = .2, fill = 'red')  + annotate("rect", xmin = 799, xmax = 921, ymin = max(conservation_CPS_I_II) - max(conservation_CPS_I_II)*0.9, ymax = max(conservation_CPS_I_II)*0.9, alpha = .2, fill = 'cyan') + annotate("rect", xmin = 1460, xmax = 1806, ymin = 0, ymax = max(conservation_CPS_I_II), alpha = .2, fill = 'purple') + annotate("rect", xmin = 1812, xmax = 1907, ymin = 0, ymax = max(conservation_CPS_I_II), alpha = .2, fill = 'yellow1') + annotate("rect", xmin = 1925, xmax = 2065, ymin = 0, ymax = max(conservation_CPS_I_II), alpha = .4, fill = 'orchid')# + annotate("rect", xmin = 1920, xmax = 2224, ymin = 0, ymax = max(conservation_CPS_I), alpha = .2, fill = 'green1') 
p5 <- ggplot() + geom_line(aes(vecI_III,conservation_CPS_I_III)) + annotate("rect", xmin = 55, xmax = 409, ymin = 0, ymax = max(conservation_CPS_I), alpha = .2, fill = 'blue') + annotate("rect", xmin = 430, xmax = 1502, ymin = 0, ymax = max(conservation_CPS_I), alpha = .2, fill = 'orange') + annotate("rect", xmin = 1371, xmax = 1486, ymin = max(conservation_CPS_I) - max(conservation_CPS_I)*0.9, ymax = max(conservation_CPS_I)*0.9, alpha = .2, fill = 'red')  + annotate("rect", xmin = 850, xmax = 973, ymin = max(conservation_CPS_I) - max(conservation_CPS_I)*0.9, ymax = max(conservation_CPS_I)*0.9, alpha = .2, fill = 'cyan')
p6 <- ggplot() + geom_line(aes(vecII_III,conservation_CPS_II_III)) + annotate("rect", xmin = 1, xmax = 355, ymin = 0, ymax = max(conservation_CPS_II), alpha = .2, fill = 'blue') + annotate("rect", xmin = 391, xmax = 1440, ymin = 0, ymax = max(conservation_CPS_II), alpha = .2, fill = 'orange') + annotate("rect", xmin = 1313, xmax = 1438, ymin = max(conservation_CPS_II) - max(conservation_CPS_II)*0.9, ymax = max(conservation_CPS_II)*0.9, alpha = .2, fill = 'red')  + annotate("rect", xmin = 799, xmax = 921, ymin = max(conservation_CPS_II) - max(conservation_CPS_II)*0.9, ymax = max(conservation_CPS_II)*0.9, alpha = .2, fill = 'cyan') + annotate("rect", xmin = 1460, xmax = 1806, ymin = 0, ymax = max(conservation_CPS_II), alpha = .2, fill = 'purple') + annotate("rect", xmin = 1812, xmax = 1907, ymin = 0, ymax = max(conservation_CPS_II), alpha = .2, fill = 'yellow1') + annotate("rect", xmin = 1925, xmax = 2065, ymin = 0, ymax = max(conservation_CPS_II), alpha = .4, fill = 'orchid')# + annotate("rect", xmin = 1920, xmax = 2224, ymin = 0, ymax = max(conservation_CPS_I), alpha = .2, fill = 'green1') 
p7 <- ggplot() + geom_line(aes(vecI_II_III,conservation_CPS_I_II_III)) + annotate("rect", xmin = 1, xmax = 355, ymin = 0, ymax = max(conservation_CPS_I_II), alpha = .2, fill = 'blue') + annotate("rect", xmin = 391, xmax = 1440, ymin = 0, ymax = max(conservation_CPS_I_II), alpha = .2, fill = 'orange') + annotate("rect", xmin = 1313, xmax = 1438, ymin = max(conservation_CPS_I_II) - max(conservation_CPS_I_II)*0.9, ymax = max(conservation_CPS_I_II)*0.9, alpha = .2, fill = 'red')  + annotate("rect", xmin = 799, xmax = 921, ymin = max(conservation_CPS_I_II) - max(conservation_CPS_I_II)*0.9, ymax = max(conservation_CPS_I_II)*0.9, alpha = .2, fill = 'cyan') + annotate("rect", xmin = 1460, xmax = 1806, ymin = 0, ymax = max(conservation_CPS_I_II), alpha = .2, fill = 'purple') + annotate("rect", xmin = 1812, xmax = 1907, ymin = 0, ymax = max(conservation_CPS_I_II), alpha = .2, fill = 'yellow1') + annotate("rect", xmin = 1925, xmax = 2065, ymin = 0, ymax = max(conservation_CPS_I_II), alpha = .4, fill = 'orchid')# + annotate("rect", xmin = 1920, xmax = 2224, ymin = 0, ymax = max(conservation_CPS_I), alpha = .2, fill = 'green1') 
p8 <- ggplot() + geom_line(aes(vec_cons, conservation_consensus)) + annotate("rect", xmin = 1, xmax = 355, ymin = 0, ymax = max(conservation_consensus), alpha = .2, fill = 'blue') + annotate("rect", xmin = 391, xmax = 1440, ymin = 0, ymax = max(conservation_consensus), alpha = .2, fill = 'orange') + annotate("rect", xmin = 1313, xmax = 1438, ymin = max(conservation_consensus) - max(conservation_consensus)*0.9, ymax = max(conservation_consensus)*0.9, alpha = .2, fill = 'red')  + annotate("rect", xmin = 799, xmax = 921, ymin = max(conservation_consensus) - max(conservation_consensus)*0.9, ymax = max(conservation_consensus)*0.9, alpha = .2, fill = 'cyan') + annotate("rect", xmin = 1460, xmax = 1806, ymin = 0, ymax = max(conservation_consensus), alpha = .2, fill = 'purple') + annotate("rect", xmin = 1812, xmax = 1907, ymin = 0, ymax = max(conservation_consensus), alpha = .2, fill = 'yellow1') + annotate("rect", xmin = 1925, xmax = 2065, ymin = 0, ymax = max(conservation_consensus), alpha = .4, fill = 'orchid')# + annotate("rect", xmin = 1920, xmax = 2224, ymin = 0, ymax = max(conservation_CPS_I), alpha = .2, fill = 'green1') 

lay <- rbind(c(1,2,3),
             c(4,5,6),
             c(7,7,7),
             c(8,8,8))

grid.arrange(grobs = list(p1, p2, p3, p4, p5, p6, p7, p8), 
             layout_matrix = lay)

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

###APES###
library(ape)

seq_I_aligned <- msaConvert(I_aligned, type="seqinr::alignment")
seq_II_aligned <- msaConvert(II_aligned, type="seqinr::alignment")
seq_III_aligned <- msaConvert(III_aligned, type="seqinr::alignment")

seq_I_II_aligned <- msaConvert(I_II_aligned, type="seqinr::alignment")
seq_I_III_aligned <- msaConvert(I_III_aligned, type="seqinr::alignment")
seq_II_III_aligned <- msaConvert(II_III_aligned, type="seqinr::alignment")
seq_I_II_III_aligned <- msaConvert(I_II_III_aligned, type="seqinr::alignment")

d_CPSI <- dist.alignment(seq_I_aligned, "identity")
d_CPSII <- dist.alignment(seq_II_aligned, "identity")
d_CPSIII <- dist.alignment(seq_III_aligned, "identity")

d_CPS_I_II <- dist.alignment(seq_I_II_aligned, "identity")
d_CPS_I_III <- dist.alignment(seq_I_III_aligned, "identity")
d_CPS_II_III <- dist.alignment(seq_II_III_aligned, "identity")
d_CPS_I_II_III <- dist.alignment(seq_I_II_III_aligned, "identity")

t_CPSI <- nj(d_CPSI)
t_CPSII <- nj(d_CPSII)
t_CPSIII <- nj(d_CPSIII)

t_CPS_I_II <- nj(d_CPS_I_II)
t_CPS_I_III <- nj(d_CPS_I_III)
t_CPS_II_III <- nj(d_CPS_II_III)
t_CPS_I_II_III <- nj(d_CPS_I_II_III)

write.tree(t_CPSI, file = "t_CPSI.nwk", append = FALSE,
           digits = 10, tree.names = TRUE)

write.tree(t_CPSII, file = "t_CPSII.nwk", append = FALSE,
           digits = 10, tree.names = TRUE)

write.tree(t_CPSIII, file = "t_CPSIII.nwk", append = FALSE,
           digits = 10, tree.names = TRUE)

write.tree(t_CPS_I_II, file = "t_CPS_I_II.nwk", append = FALSE,
           digits = 10, tree.names = TRUE)

write.tree(t_CPS_I_III, file = "t_CPS_I_III.nwk", append = FALSE,
           digits = 10, tree.names = TRUE)

write.tree(t_CPS_II_III, file = "t_CPS_II_III.nwk", append = FALSE,
           digits = 10, tree.names = TRUE)

write.tree(t_CPS_I_II_III, file = "t_CPS_I_II_III.nwk", append = FALSE,
           digits = 10, tree.names = TRUE)

###APES###


msaPrettyPrint(msaData, y=c(100, 213), output="pdf",
               showNames="none", showLogo="top", askForOverwrite=FALSE)
