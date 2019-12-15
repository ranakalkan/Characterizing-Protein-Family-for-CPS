library(msa)
setwd("/home/fatih/Documents/ENS210/project/fastas/specific")

amph_1 <- readAAStringSet("CPS1_amphibians.fas")
fish_1 <- readAAStringSet("CPS1_fishes.fas")
ter_1 <- readAAStringSet("CPS1_birds&mammals.fas")

amph_2 <- readAAStringSet("CPS2_amphibians.fas")
fish_2 <- readAAStringSet("CPS2_fishes.fas")
ter_2 <- readAAStringSet("CPS2_terrastrial.fas")

amph_1_a <- msa(amph_1)
fish_1_a <- msa(fish_1)
ter_1_a <- msa(ter_1)

amph_2_a <- msa(amph_2)
fish_2_a <- msa(fish_2)
ter_2_a <- msa(ter_2)

writeXStringSet(unmasked(amph_1_a), filepath = "amph_1_a.fas")
writeXStringSet(unmasked(fish_1_a), filepath = "fish_1_a.fas")
writeXStringSet(unmasked(ter_1_a), filepath = "ter_1_a.fas")

writeXStringSet(unmasked(amph_2_a), filepath = "amph_2_a.fas")
writeXStringSet(unmasked(fish_2_a), filepath = "fish_2_a.fas")
writeXStringSet(unmasked(ter_2_a), filepath = "ter_2_a.fas")

amph_1_c <- readAAStringSet("amph_1_a_c.fas")
fish_1_c <- readAAStringSet("fish_1_a_c.fas")
ter_1_c <- readAAStringSet("ter_1_a_c.fas")

amph_2_c <- readAAStringSet("amph_2_a_c.fas")
fish_2_c <- readAAStringSet("fish_2_a_c.fas")
ter_2_c <- readAAStringSet("ter_2_a_c.fas")

amph_1_a_c <- msa(amph_1_c)
fish_1_a_c <- msa(fish_1_c)
ter_1_a_c <- msa(ter_1_c)

amph_2_a_c <- msa(amph_2_c)
fish_2_a_c <- msa(fish_2_c)
ter_2_a_c <- msa(ter_2_c)

data(BLOSUM62)

conservation_fish_1 = msa::msaConservationScore(fish_1_a_c, BLOSUM62)
conservation_amph_1 = msa::msaConservationScore(amph_1_a_c, BLOSUM62)
conservation_ter_1 = msa::msaConservationScore(ter_1_a_c, BLOSUM62)

conservation_fish_2 = msa::msaConservationScore(fish_2_a_c, BLOSUM62)
conservation_amph_2 = msa::msaConservationScore(amph_2_a_c, BLOSUM62)
conservation_ter_2 = msa::msaConservationScore(ter_2_a_c, BLOSUM62)

pos1 <- vector()
pos1 <- c(pos1, 1:length(conservation_fish_1))

pos2 <- vector()
pos2 <- c(pos2, 1:length(conservation_amph_1))

pos3 <- vector()
pos3 <- c(pos3, 1:length(conservation_ter_1))

pos4 <- vector()
pos4 <- c(pos4, 1:length(conservation_fish_2))

pos5 <- vector()
pos5 <- c(pos5, 1:length(conservation_amph_2))

pos6 <- vector()
pos6 <- c(pos6, 1:length(conservation_ter_2))


p1 <- ggplot() + geom_line(aes(pos1,conservation_fish_1)) + annotate("rect", xmin = 55, xmax = 409, ymin = 0, ymax = max(conservation_fish_1), alpha = .2, fill = 'blue') + annotate("rect", xmin = 430, xmax = 1502, ymin = 0, ymax = max(conservation_fish_1), alpha = .2, fill = 'orange') + annotate("rect", xmin = 1371, xmax = 1486, ymin = max(conservation_fish_1) - max(conservation_fish_1)*0.9, ymax = max(conservation_fish_1)*0.9, alpha = .2, fill = 'red')  + annotate("rect", xmin = 850, xmax = 973, ymin = max(conservation_fish_1) - max(conservation_fish_1)*0.9, ymax = max(conservation_fish_1)*0.9, alpha = .2, fill = 'cyan')
p2 <- ggplot() + geom_line(aes(pos2,conservation_amph_1)) + annotate("rect", xmin = 55, xmax = 409, ymin = 0, ymax = max(conservation_amph_1), alpha = .2, fill = 'blue') + annotate("rect", xmin = 430, xmax = 1502, ymin = 0, ymax = max(conservation_amph_1), alpha = .2, fill = 'orange') + annotate("rect", xmin = 1371, xmax = 1486, ymin = max(conservation_amph_1) - max(conservation_amph_1)*0.9, ymax = max(conservation_amph_1)*0.9, alpha = .2, fill = 'red')  + annotate("rect", xmin = 850, xmax = 973, ymin = max(conservation_amph_1) - max(conservation_amph_1)*0.9, ymax = max(conservation_amph_1)*0.9, alpha = .2, fill = 'cyan')
p3 <- ggplot() + geom_line(aes(pos3,conservation_ter_1)) + annotate("rect", xmin = 55, xmax = 409, ymin = 0, ymax = max(conservation_ter_1), alpha = .2, fill = 'blue') + annotate("rect", xmin = 430, xmax = 1502, ymin = 0, ymax = max(conservation_ter_1), alpha = .2, fill = 'orange') + annotate("rect", xmin = 1371, xmax = 1486, ymin = max(conservation_ter_1) - max(conservation_ter_1)*0.9, ymax = max(conservation_ter_1)*0.9, alpha = .2, fill = 'red')  + annotate("rect", xmin = 850, xmax = 973, ymin = max(conservation_ter_1) - max(conservation_ter_1)*0.9, ymax = max(conservation_ter_1)*0.9, alpha = .2, fill = 'cyan')

p4 <- ggplot() + geom_line(aes(pos4,conservation_fish_2)) + annotate("rect", xmin = 1, xmax = 355, ymin = 0, ymax = max(conservation_fish_2), alpha = .2, fill = 'blue') + annotate("rect", xmin = 391, xmax = 1440, ymin = 0, ymax = max(conservation_fish_2), alpha = .2, fill = 'orange') + annotate("rect", xmin = 1313, xmax = 1438, ymin = max(conservation_fish_2) - max(conservation_fish_2)*0.9, ymax = max(conservation_fish_2)*0.9, alpha = .2, fill = 'red')  + annotate("rect", xmin = 799, xmax = 921, ymin = max(conservation_fish_2) - max(conservation_fish_2)*0.9, ymax = max(conservation_fish_2)*0.9, alpha = .2, fill = 'cyan') + annotate("rect", xmin = 1460, xmax = 1806, ymin = 0, ymax = max(conservation_fish_2), alpha = .2, fill = 'purple') + annotate("rect", xmin = 1812, xmax = 1907, ymin = 0, ymax = max(conservation_fish_2), alpha = .2, fill = 'yellow1') + annotate("rect", xmin = 1925, xmax = 2065, ymin = 0, ymax = max(conservation_fish_2), alpha = .4, fill = 'orchid')# + annotate("rect", xmin = 1920, xmax = 2224, ymin = 0, ymax = max(conservation_CPS_I), alpha = .2, fill = 'green1')
p5 <- ggplot() + geom_line(aes(pos5,conservation_amph_2)) + annotate("rect", xmin = 1, xmax = 355, ymin = 0, ymax = max(conservation_amph_2), alpha = .2, fill = 'blue') + annotate("rect", xmin = 391, xmax = 1440, ymin = 0, ymax = max(conservation_amph_2), alpha = .2, fill = 'orange') + annotate("rect", xmin = 1313, xmax = 1438, ymin = max(conservation_amph_2) - max(conservation_amph_2)*0.9, ymax = max(conservation_amph_2)*0.9, alpha = .2, fill = 'red')  + annotate("rect", xmin = 799, xmax = 921, ymin = max(conservation_amph_2) - max(conservation_amph_2)*0.9, ymax = max(conservation_amph_2)*0.9, alpha = .2, fill = 'cyan') + annotate("rect", xmin = 1460, xmax = 1806, ymin = 0, ymax = max(conservation_amph_2), alpha = .2, fill = 'purple') + annotate("rect", xmin = 1812, xmax = 1907, ymin = 0, ymax = max(conservation_amph_2), alpha = .2, fill = 'yellow1') + annotate("rect", xmin = 1925, xmax = 2065, ymin = 0, ymax = max(conservation_amph_2), alpha = .4, fill = 'orchid')# + annotate("rect", xmin = 1920, xmax = 2224, ymin = 0, ymax = max(conservation_CPS_I), alpha = .2, fill = 'green1')
p6 <- ggplot() + geom_line(aes(pos6,conservation_ter_2)) + annotate("rect", xmin = 1, xmax = 355, ymin = 0, ymax = max(conservation_ter_2), alpha = .2, fill = 'blue') + annotate("rect", xmin = 391, xmax = 1440, ymin = 0, ymax = max(conservation_ter_2), alpha = .2, fill = 'orange') + annotate("rect", xmin = 1313, xmax = 1438, ymin = max(conservation_ter_2) - max(conservation_ter_2)*0.9, ymax = max(conservation_ter_2)*0.9, alpha = .2, fill = 'red')  + annotate("rect", xmin = 799, xmax = 921, ymin = max(conservation_ter_2) - max(conservation_ter_2)*0.9, ymax = max(conservation_ter_2)*0.9, alpha = .2, fill = 'cyan') + annotate("rect", xmin = 1460, xmax = 1806, ymin = 0, ymax = max(conservation_ter_2), alpha = .2, fill = 'purple') + annotate("rect", xmin = 1812, xmax = 1907, ymin = 0, ymax = max(conservation_ter_2), alpha = .2, fill = 'yellow1') + annotate("rect", xmin = 1925, xmax = 2065, ymin = 0, ymax = max(conservation_ter_2), alpha = .4, fill = 'orchid')# + annotate("rect", xmin = 1920, xmax = 2224, ymin = 0, ymax = max(conservation_CPS_I), alpha = .2, fill = 'green1')

lay <- rbind(c(1,1,1),
             c(2,2,2),
             c(3,3,3),
             c(4,4,4),
             c(5,5,5),
             c(6,6,6))

grid.arrange(grobs = list(p1, p2, p3, p4, p5, p6), 
             layout_matrix = lay)




