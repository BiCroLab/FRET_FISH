setwd("/Users/anamota/Documents/projects/FRET-FISH/revision/Lamina dist and 3D segmentation/LADs")


library(tidyr)
library(ggplot2)
library(dplyr)
library(ggpubr)

data_wide = read.delim("FRET_eff.txt", header = T, as.is = T, sep = ",")



my_comparisons <- list( c("Minar2.cLAD", "Grxcr2.cLAD"), c("Minar2.cLAD", "X4930426D05Rik.cLAD"), 
                        c("Minar2.cLAD", "Hspa9.iLAD"), c("Minar2.cLAD", "Nars.iLAD"),
                        c("Minar2.cLAD", "Atp5a1.iLAD"), c("Minar2.cLAD", "Pkm.chr9"),
                        c("Minar2.cLAD", "P4hb.chr11"), c("Minar2.cLAD", "Hsp90ab1.chr17"))


countX <- NULL
countX$Gene = c("Minar2.cLAD","Grxcr2.cLAD","X4930426D05Rik.cLAD","Hspa9.iLAD","Nars.iLAD","Atp5a1.iLAD", "Pkm.chr9", "P4hb.chr11", "Hsp90ab1.chr17")
countX$ypos = c(65, 65, 65, 65, 65, 65, 65, 65, 65)
countX$Gene <- factor(countX$Gene, levels = c("Minar2.cLAD","Grxcr2.cLAD","X4930426D05Rik.cLAD","Hspa9.iLAD","Nars.iLAD","Atp5a1.iLAD", "Pkm.chr9", "P4hb.chr11", "Hsp90ab1.chr17"))
countX = as.data.frame(countX)


setwd("/Users/anamota/Desktop")
pdf("FRET.pdf", height=4, width=6) 

data_wide %>% 
  gather(key="Gene", value="FRET_efficiency") %>%
  mutate( Gene=factor(Gene,levels=c("Minar2.cLAD", "Grxcr2.cLAD", "X4930426D05Rik.cLAD", "Hspa9.iLAD", "Nars.iLAD", "Atp5a1.iLAD", "Pkm.chr9", "P4hb.chr11", "Hsp90ab1.chr17")) ) %>%
  ggplot( aes(x=Gene, y=FRET_efficiency, color=factor(Gene))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  ylab("FRET Efficiency (%)") + xlab("") +
  scale_color_manual(values=c("#E69F00","#E69F00","#E69F00","#56B4E9","#56B4E9","#56B4E9", "#B22222","#B22222","#B22222")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text=element_text(family="Helvetica", size=7)) +
  geom_text(data = countX, mapping = aes(x=Gene, y=ypos, label=paste("n = ", colSums(!is.na(data_wide)))), size=2, colour="black") +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", family="Helvetica", size=2) # Add pairwise comparisons p-value


dev.off()