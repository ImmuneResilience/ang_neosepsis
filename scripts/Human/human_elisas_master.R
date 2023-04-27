######################
# Vascular resilience
# Ang1 and Ang2: master human ELISA analysis
#####################

# In this script, we will generate all ELISA boxplots generated for the manuscripts for human data, Angiopoetin 1 and Angiopoetin 2.

# load data
elisa.master <- read.csv("Rdata/Human/human_elisa.csv")

# load packages
library(plyr)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)

# set output
output.dir <- "Figures/human/boxplots/"

# Set plotting aesthetics -------------------------------------------------
box.theme <- theme(panel.grid = element_blank(),
                   axis.text = element_text(size = 13, colour = "black"),
                   strip.text = element_text(size = 14),
                   axis.title = element_text(size = 14),
                   axis.title.x = element_blank(),
                   legend.position = "none",
                   plot.title = element_text(size = 14, hjust = 0.5))

cor.theme <- theme(panel.grid = element_blank(),
                   axis.text = element_text(size = 11, colour = "black"),
                   strip.text = element_text(size = 14),
                   axis.title = element_text(size = 14),
                   legend.position = "none",
                   axis.text.x = element_text(size = 13),
                   plot.title = element_text(size = 14, hjust = 0.5))

# Set Colours
cols <- brewer.pal("Dark2", n = 3)
cols <- cols[c(2:3)]



# Generating Boxplots -----------------------------------------------


# Ang1_Ang2_elisas_human_clinical: boxplot ####
elisa_ang1_ang2 <- ggplot(elisa.master, aes(group, ng_ml, fill = group)) +
  geom_violin(draw_quantiles = TRUE, alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, size = 0.75, width = 0.3) +
  geom_point(position = position_jitter(width = 0.05), shape = 21, size = 0.8) +
  stat_compare_means(size = 4, method = "wilcox", comparisons = list(c("Controls", "Cases")), label = "p.signif") +
  labs(y = "Plasma concentration (ng per ml)") +
  facet_wrap(~protein, scales = "free") +
  scale_fill_manual(values = cols) +
  theme_classic() + 
  box.theme

elisa_ang1_ang2

# print
ggsave(paste0(output.dir, "violin_human_elisa_case_control_a1a2.pdf"),
       device = "pdf",
       dpi = 300,
       width = 4, height = 3.5)

#sPLA2_elisa_controls_vs_cases ####

elisa_sPLA2 <- ggplot(elisa.master, aes(group, log2_sPLA2, fill = group)) +
  geom_violin(draw_quantiles = TRUE, alpha = 0.5, width = 1.1) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, size = 0.75, width = 0.18) +
  geom_point(position = position_jitter(width = 0.05), shape = 21, size = 0.8) +
  stat_compare_means(method = "wilcox", label = "p.signif", comparisons = list(c("Controls", "Cases"))) +
  labs(y = "Concentration (log2 pg/ml)", title = "Plasma sPLA2") + 
  scale_fill_manual(values = cols) +
  theme_bw() + 
  box.theme +
  ylim(c(10,22))

elisa_sPLA2

# print
ggsave(paste0(output.dir, "violin_human_elisa_case_control_sPLA2.pdf"),
       device = "pdf",
       dpi = 300,
       width = 2.5, height = 3.5)

  
# ratio_Ang1_Ang2_elisas_human_clinical ####

elisa_ang1ang2_ratio <- ggplot(elisa.master, aes(group, ratio, fill = group)) +
  geom_violin(draw_quantiles = TRUE, alpha = 0.5, width = 0.7) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, size = 0.75, width = 0.18) +
  geom_point(position = position_jitter(width = 0.05), shape = 21, size = 0.8) +
  stat_compare_means(method = "wilcox", label = "p.signif", comparisons = list(c("Controls", "Cases"))) +
  labs(title = "Ang1:Ang2 Ratio", y = "Ratio") + 
  scale_fill_manual(values = cols) +
  theme_bw() + 
  box.theme +
  ylim(c(0,4.5))

elisa_ang1ang2_ratio

# print
ggsave(paste0(output.dir, "violin_human_elisa_case_control_a1a2_ratio.pdf"),
       device = "pdf",
       dpi = 300,
       width = 2.5, height = 3.5)


# Correlation plots -------------------------------------------------------


# Ang1 vs. sPLA2 ####
a1 <- elisa.master[elisa.master$protein == "Ang1",]

test <- cor.test(~ng_ml +  log2_sPLA2, method = "pearson", data = a1)

# get results
pval <- round(test$p.value, digits = 4)
r <- round(test$estimate, digits = 3)

# create label
res.lab <- paste0("Pearson's r = ", r, ", p = ", pval)
res.lab

# plot

ggplot(a1, aes(x = ng_ml, y = log2_sPLA2)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  cor.theme +
  labs(x = "Plasma Ang1 (ng/mL)", y = "Plasma sPLA2 (Log2 pg/mL)", title = "Ang1 vs sPLA2") +
  annotate("text",label = res.lab, x = 25, y = 25)

# save
ggsave(paste0(output.dir, "human_corplot_a1vspla2.pdf"),
       device = "pdf",
       dpi = 300,
       width = 3.5, height = 3.5)


# Ang2 vs. sPLA2 ####
a2 <- elisa.master[elisa.master$protein == "Ang2",]

test <- cor.test(~ng_ml +  log2_sPLA2, method = "pearson", data = a2)

# get results
pval <- round(test$p.value, digits = 4)
r <- round(test$estimate, digits = 3)

# create label
res.lab <- paste0("Pearson's r = ", r, ", p = ", pval)
res.lab

# plot

ggplot(a2, aes(x = ng_ml, y = log2_sPLA2)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  cor.theme +
  labs(x = "Plasma Ang2 (ng/mL)", y = "Plasma sPLA2 (Log2 pg/mL)", title = "Ang2 vs sPLA2") +
  annotate("text",label = res.lab, x = 12, y = 23)

# save
ggsave(paste0(output.dir, "human_corplot_a2vspla2.pdf"),
       device = "pdf",
       dpi = 300,
       width = 3.5, height = 3.5)

# Ang1/2 ratio vs. sPLA2 ####

test <- cor.test(~ ratio +  log2_sPLA2, method = "pearson", data = a1)

# get results
pval <- round(test$p.value, digits = 4)
r <- round(test$estimate, digits = 3)

# create label
res.lab <- paste0("Pearson's r = ", r, ", p = ", pval)
res.lab

# plot

ggplot(a1, aes(x = ratio, y = log2_sPLA2)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  cor.theme +
  labs(x = "Plasma Ang1:2 ratio", y = "Plasma sPLA2 (Log2 pg/mL)", title = "Ang1:2 ratio vs sPLA2") +
  annotate("text",label = res.lab, x = 2, y = 23)

# save
ggsave(paste0(output.dir, "human_corplot_a1a2ratiovspla2.pdf"),
       device = "pdf",
       dpi = 300,
       width = 3.5, height = 3.5)

