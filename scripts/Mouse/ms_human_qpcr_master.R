######################
# Vascular resilience
# Human qPCR: master analysis
#####################

# In this script, we will generate all qPCR boxplots and correlations for human data.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# install packages
library(tidyverse)
library(ggpubr)
library(RColorBrewer)

# load data
dat <- read.csv("../../Rdata/Human/ms_human_qpcr_master.csv")

# set output
output.dir <- "../../Figures/human/boxplots/"

# Set Theme ---------------------------------------------------------------

box.theme <- theme(panel.grid = element_blank(),
                   axis.text = element_text(size = 13, colour = "black"),
                   strip.text = element_text(size = 14),
                   axis.title = element_text(size = 14),
                   axis.title.x = element_blank(),
                   legend.position = "none",
                   axis.text.x = element_text(size = 13),
                   plot.title = element_text(size = 14, hjust = 0.5))

# Set Colours
cols <- brewer.pal("Dark2", n = 3)
cols <- cols[c(2:3)]


# Analyze data ------------------------------------------------------------

# Set to Sepsis/Control, to be consistent with ELISA data
dat$Type <- factor(dat$Type, levels = c("Sepsis", "Control"))

# Plot Alox5, CYP2J2, Ptgs ####

ggplot(dat, aes(Type, logfc, fill = Type)) +
  #geom_violin(draw_quantiles = TRUE, alpha = 0.5, width = 2) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, size = 0.75, width = 0.6) +
  geom_point(position = position_jitter(width = 0.05), shape = 21, size = 1.5) +
  stat_compare_means(size = 3, comparisons = list(c("Control", "Sepsis")), label = "p.format") +
  facet_wrap(~gene, scales = "free", ncol = 4) +
  scale_fill_manual(values = cols) +
  labs(y = "Log Fold-Change") +
  theme_classic() + 
  box.theme


# save
ggsave(paste0(output.dir, "human_qpcr_control_case_aloxCYPPtgs.pdf"),
       device = "pdf",
       dpi = 300,
       width = 5.4, height = 3.5)

# END ####