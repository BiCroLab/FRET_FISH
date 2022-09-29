## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for 

## Load/install packages
packages = c("data.table", "rtracklayer", "ggrepel")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

# Load in data
dt = import.bw("/mnt/sequencing/FRET-FISH_HiC/atacseq/data/GSM3656632_7F-D0-rep1.bw") # Publicly available
# dt = import.bw("/mnt/sequencing/FRET-FISH_HiC/atacseq/nextfdlow/bwa/mergedLibrary/bigwig/one_R1.mLb.clN.bigWig") # Own processed
probes = fread("/mnt/sequencing/FRET-FISH_HiC/FFS/probes_all.tsv")

dt = as.data.table(dt)

# Set keys
setkey(dt, seqnames, start, end)
setkey(probes, chr, start, end)

# Overlap
overlap = foverlaps(probes, dt)

# sum_score = overlap[, .(sum_atacseq = mean(score_norm)), by = .(name, fret_eff)]
sum_score = overlap[, .(mean_atacseq = mean(score)), by = .(name, fret_eff)]

sum_score_x = sum_score[mean_score$name %in% probes[chr == "chrX", name]]
sum_score_chr18 = sum_score[mean_score$name %in% probes[chr == "chr18", name]]

# Plot
plt1 = ggplot(sum_score, aes(x = fret_eff, y = mean_atacseq, label = name)) +
  geom_smooth(method='lm', se = F, color = "red", size = 1, linetype = 2) +
  geom_point(size = 4, color = "lightblue") +
  geom_text_repel() +
  labs(y = "Mean ATAC-seq score", x = "FFS") +
  stat_cor(method = "pearson", size = 4, aes(label = paste0("PCC: ", ..label..))) +
  stat_cor(method = "spearman", size = 4, position = position_nudge(y = 1), aes(label = paste0("SCC: ", ..label..)))

plt2 = ggplot(sum_score_x, aes(x = fret_eff, y = mean_atacseq, label = name)) +
  geom_smooth(method='lm', se = F, color = "red", size = 1, linetype = 2) +
  geom_point(size = 4, color = "lightblue") +
  geom_text_repel() +
  labs(y = "Mean ATAC-seq score", x = "FFS") +
  stat_cor(method = "pearson", size = 4, aes(label = paste0("PCC: ", ..label..))) +
  stat_cor(method = "spearman", size = 4, position = position_nudge(y = 500), aes(label = paste0("SCC: ", ..label..)))

plt3 = ggplot(sum_score_chr18, aes(x = fret_eff, y = mean_atacseq, label = name)) +
  geom_smooth(method='lm', se = F, color = "red", size = 1, linetype = 2) +
  geom_point(size = 4, color = "lightblue") +
  geom_text_repel() +
  labs(y = "Mean ATAC-seq score", x = "FFS") +
  stat_cor(method = "pearson", size = 4, aes(label = paste0("PCC: ", ..label..))) +
  stat_cor(method = "spearman", size = 4, position = position_nudge(y = 500), aes(label = paste0("SCC: ", ..label..)))

# # Save
# save_and_plot(plt1, "/mnt/sequencing/FRET-FISH_HiC/plots/ATACseq/ATAC_FFS_all", pointsize = 1, width = 7, height = 7)
# save_and_plot(plt2, "/mnt/sequencing/FRET-FISH_HiC/plots/ATACseq/ATAC_FFS_chrx", pointsize = 1, width = 7, height = 7)
# save_and_plot(plt3, "/mnt/sequencing/FRET-FISH_HiC/plots/ATACseq/ATAC_FFS_chr18", pointsize = 1, width = 7, height = 7)
