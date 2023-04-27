######################
# Vascular resilience
# Ang1 and Ang2: master ELISA analysis
#####################

# In this script, we will generate all ELISA boxplots generated for the manuscripts for mouse data, Angiopoetin 1 and Angiopoetin 2.


# Setup -------------------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load data
elisa.master <- read.csv("../../Rdata/Mouse/mouse_elisa_masterfile.csv")

# load packages
library(plyr)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)

# set output
output.dir <- "../../Figures/mouse/boxplots/"

# Set plotting aesthetics -------------------------------------------------

box.theme <- theme(panel.grid = element_blank(),
                   axis.text = element_text(size = 11, colour = "black"),
                   strip.text = element_text(size = 14),
                   axis.title = element_text(size = 14),
                   axis.title.x = element_blank(),
                   legend.position = "none",
                   axis.text.x = element_text(size = 13),
                   plot.title = element_text(size = 14, hjust = 0.5))


# colours
cols <- brewer.pal(n = 4, "Dark2")[2:4]
cols <- c("#969696", cols)

# Prepare variables -------------------------------------------------------


elisa.master$treatment <- gsub("cont", "control", elisa.master$treatment)

# Subset data to experiment -----------------------------------------------

elisa.list <- split(elisa.master, elisa.master$experiment)

# ARA vs CS vs Control ----------------------------------------------------

# for this analysis, we gather the aa, cs, and control groups.  Both the 2 and 4 hour post challenge groups are used.

aa <- elisa.list$aa.lname.arg_serum %>%
  filter(treatment %in% c("aa", "cs", "control"))

# set variable names

aa$treatment <- gsub("aa", "aa+cs", aa$treatment)

# set statistical comparisons
aa_comparisons <-  list(c("control", "cs"),
                        c("cs", "aa+cs"),
                        c("control", "aa+cs"))

# create fold change from controls
aa <- ddply(aa,
            .(analyte),
            transform,
            fold_change = concentration_ng.mL/median(concentration_ng.mL[treatment == "control"]) - 1)


# set factor levels
aa$treatment <- factor(aa$treatment, levels = c("control",
                                                "cs",
                                                "aa+cs"))

# plot concentration
ggplot(aa, aes(treatment, concentration_ng.mL, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, size = 1) +
  geom_point(position = position_jitter(0.15), shape = 21, size = 2, alpha = 0.6) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = aa_comparisons, size = 4) +
  facet_wrap(~analyte, scales = "free") +
  scale_fill_manual(values = cols) +
  #ggtitle("AA Serum Ang1 and Ang2") +
  theme_classic() + 
  box.theme +
  ylab("Concentration (ng/mL)") +
  scale_color_manual(values = cols)

# save plot
ggsave(paste0(output.dir, "mouse_elisa_aa_cs_control_a1a2.pdf"),
              device = "pdf",
              dpi = 300,
              width = 4.5, height = 3.5)


# plot fold change
ggplot(aa, aes(treatment, fold_change, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, size = 1) +
  geom_point(position = position_jitter(0.15), shape = 21, size = 2, alpha = 0.6) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = aa_comparisons, size = 4) +
  facet_wrap(~analyte, scales = "free") +
  scale_fill_manual(values = cols) +
  #ggtitle("AA Serum Ang1 and Ang2") +
  theme_classic() + 
  box.theme +
  ylab("Fold Change") +
  scale_color_manual(values = cols)

ggsave(paste0(output.dir, "mouse_elisa_aa_cs_control_a1a2_FC.pdf"),
       device = "pdf",
       dpi = 300,
       width = 4.5, height = 3.5)

# LNAME vs. CS vs. Control -------------------------------------------------------
# for this analysis, we gather the lname, cs, and control groups.  Only the 2 hour post challenge groups are used.

ln <- elisa.list$aa.lname.arg_serum %>%
  filter(treatment %in% c("lname", "cs", "control") | is.na(hpc),
         is.na(hpc) | hpc == 2)

table(ln$treatment)

# set variable names
ln$treatment <- gsub("lname", "lname+cs", ln$treatment)

# create fold change from controls
ln <- ddply(ln,
            .(analyte),
            transform,
            fold_change = concentration_ng.mL/median(concentration_ng.mL[treatment == "control"]) - 1)

# set comparions
ln_comparisons <-  list(c("control", "cs"),
                        c("cs", "lname+cs"),
                        c("control", "lname+cs"))

# set factor levels
ln$treatment <- factor(ln$treatment, levels = c("control", "cs", "lname+cs"))

# plot concentration
ggplot(ln, aes(treatment, concentration_ng.mL, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, size = 1) +
  geom_point(position = position_jitter(0.15), shape = 21, size = 2, alpha = 0.6) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = ln_comparisons, size = 4) +
  facet_wrap(~analyte, scales = "free") +
  scale_fill_manual(values = cols) +
  #ggtitle("LNAME Serum Ang1 and Ang2") +
  theme_classic() + 
  box.theme +
  ylab("Concentration (ng/mL)")

# save plot
ggsave(paste0(output.dir, "mouse_elisa_lname_cs_control_a1a2.pdf"),
       device = "pdf",
       dpi = 300,
       width = 4.5, height = 3.5)

# plot fold change
ggplot(ln, aes(treatment, fold_change, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, size = 1) +
  geom_point(position = position_jitter(0.15), shape = 21, size = 2, alpha = 0.6) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = ln_comparisons, size = 4) +
  facet_wrap(~analyte, scales = "free") +
  scale_fill_manual(values = cols) +
  #ggtitle("LNAME Serum Ang1 and Ang2") +
  theme_classic() + 
  box.theme +
  ylab("Fold Change")

# save plot
ggsave(paste0(output.dir, "mouse_elisa_lname_cs_control_a1a2_FC.pdf"),
       device = "pdf",
       dpi = 300,
       width = 4.5, height = 3.5)

# Ang1 + L-Arg ------------------------------------------------------------

# for this experiment, we gather the a1+Larg experiment.
a1la <- elisa.list$`a1+arg_serum`

# create fold change from controls
a1la <- ddply(a1la,
            .(analyte),
            transform,
            fold_change = concentration_ng.mL/median(concentration_ng.mL[treatment == "control"]) - 1)

# set comparions
a1la_comparisons <-  list(c("control", "cs"),
                        c("control", "a1+arg+cs"),
                        c("control", "a1+arg"),
                        c("cs", "a1+arg+cs"),
                        c("cs", "a1+arg"),
                        c("a1+arg+cs", "a1+arg"))

# set factor levels
a1la$treatment <- factor(a1la$treatment, levels = c("control", "cs", "a1+arg+cs", "a1+arg"))

# plot concentration
ggplot(a1la, aes(treatment, concentration_ng.mL, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, size = 1) +
  geom_point(position = position_jitter(0.15), shape = 21, size = 2, alpha = 0.6) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = a1la_comparisons, size = 4) +
  facet_wrap(~analyte, scales = "free") +
  scale_fill_manual(values = cols) +
  #ggtitle("Ang1+Arg Serum Ang1") +
  theme_classic() + 
  box.theme +
  ylab("Concentration (ng/mL)") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1.0, vjust = 1.0)) +
  ylim(c(0, 420))

# save plot
ggsave(paste0(output.dir, "mouse_elisa_a1arg_cs_control_a1.pdf"),
       device = "pdf",
       dpi = 300,
       width = 3, height = 3.5)

# plot fold change
ggplot(a1la, aes(treatment, fold_change, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, size = 1) +
  geom_point(position = position_jitter(0.15), shape = 21, size = 2, alpha = 0.6) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = a1la_comparisons, size = 4) +
  facet_wrap(~analyte, scales = "free") +
  scale_fill_manual(values = cols) +
  #ggtitle("Ang1+Arg Serum Ang1") +
  theme_classic() + 
  box.theme +
  ylab("Fold Change") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1.0, vjust = 1.0))+
  ylim(c(-1, 1.6))

# save plot
ggsave(paste0(output.dir, "mouse_elisa_a1arg_cs_control_a1_FC.pdf"),
       device = "pdf",
       dpi = 300,
       width = 3, height = 3.5)


# Ang1 + LNAME vs control -------------------------------------------------

# for this experiment, we gather the a1+Larg experiment.
a1ln <- elisa.list$`a1+lname_serum`

# create fold change from controls
a1ln <- ddply(a1ln,
              .(analyte),
              transform,
              fold_change = concentration_ng.mL/median(concentration_ng.mL[treatment == "control"]) - 1)

# set comparions
a1ln_comparisons <-  list(c("control", "cs"),
                          c("control", "a1+lname+cs"),
                          c("control", "a1+lname"),
                          c("cs", "a1+lname+cs"),
                          c("cs", "a1+lname"),
                          c("a1+lname+cs", "a1+lname"))

# set factor levels
a1ln$treatment <- factor(a1ln$treatment, levels = c("control", "cs", "a1+lname+cs", "a1+lname"))

# plot concentration
ggplot(a1ln, aes(treatment, concentration_ng.mL, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, size = 1) +
  geom_point(position = position_jitter(0.15), shape = 21, size = 2, alpha = 0.6) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = a1ln_comparisons, size = 4) +
  facet_wrap(~analyte, scales = "free") +
  scale_fill_manual(values = cols) +
  #ggtitle("Ang1+LNAME Serum Ang1") +
  theme_classic() + 
  box.theme +
  ylab("Concentration (ng/mL)") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1.0, vjust = 1.0)) +
  ylim(c(0, 830))

# save plot
ggsave(paste0(output.dir, "mouse_elisa_a1lname_cs_control_a1.pdf"),
       device = "pdf",
       dpi = 300,
       width = 3, height = 3.5)

# plot fold change
ggplot(a1ln, aes(treatment, fold_change, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, size = 1) +
  geom_point(position = position_jitter(0.15), shape = 21, size = 2, alpha = 0.6) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = a1ln_comparisons, size = 4) +
  facet_wrap(~analyte, scales = "free") +
  scale_fill_manual(values = cols) +
  #ggtitle("Ang1+LNAME Serum Ang1") +
  theme_classic() + 
  box.theme +
  ylab("Fold Change") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1.0, vjust = 1.0)) +
  ylim(c(-1, 1.7))

# save plot
ggsave(paste0(output.dir, "mouse_elisa_a1lname_cs_control_a1_FC.pdf"),
       device = "pdf",
       dpi = 300,
       width = 3, height = 3.5)

