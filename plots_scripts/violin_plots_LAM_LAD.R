setwd("/Users/anamota/Documents/projects/FRET-FISH/revision/Lamina dist and 3D segmentation/LADs")


library(tidyr)
library(ggplot2)
library(dplyr)
library(ggpubr)

data_wide = read.delim("LAM_LAD.txt", header = T, as.is = T, sep = ",")



my_comparisons <- list( c("cLADs", "iLADs"))


countX <- NULL
countX$Gene = c("cLADs","iLADs")
countX$ypos = c(0, 0)
countX$Gene <- factor(countX$Gene, levels = c("cLADs","iLADs"))
countX = as.data.frame(countX)


setwd("/Users/anamota/Desktop")
pdf("LAM_LAD.pdf", height=4, width=6) 

data_wide %>% 
  gather(key="Gene", value="FRET_efficiency") %>%
  mutate( Gene=factor(Gene,levels=c("cLADs","iLADs")) ) %>%
  ggplot( aes(x=Gene, y=FRET_efficiency, color=factor(Gene))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5) +
  #geom_point(position = position_jitter(seed = 0.1, width = 0.1), size=0.2) +
  ylab("Lamina distance to center norm.") + xlab("") +
  scale_color_manual(values=c("#E69F00","#56B4E9")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text=element_text(family="Helvetica", size=7)) +
  geom_text(data = countX, mapping = aes(x=Gene, y=ypos, label=paste("n = ", colSums(!is.na(data_wide)))), size=2, colour="black") +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", family="Helvetica", size=2) # Add pairwise comparisons p-value


dev.off()