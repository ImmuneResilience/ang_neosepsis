######################
# Vascular resilience
# Ang1 and Ang2: master ELISA analysis
#####################

# In this script, we will generate all ELISA boxplots generated for the manuscripts for mouse data, Angiopoetin 1 and Angiopoetin 2.


# Setup -------------------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load data
elisa.master <- read.csv("../../Rdata/Mouse/ms_mouse_elisa_masterfile.csv")

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

# ARA vs CS vs Control ----------------------------------------------------

# for this analysis, we gather the aa, cs, and control groups.  Both the 2 and 4 hour post challenge groups are used.

aa <- filter(elisa.master, treatment %in% c("aa", "cs", "control"))

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


# plot fold change
ggplot(aa, aes(treatment, fold_change, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, size = 1) +
  geom_point(position = position_jitter(0.15), shape = 21, size = 2, alpha = 0.6) +
  stat_compare_means(method = "wilcox.test", label = "p.format", comparisons = aa_comparisons, size = 4) +
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

ln <- filter(elisa.master, treatment %in% c("lname", "cs", "control") | is.na(hpc),
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

# plot fold change
ggplot(ln, aes(treatment, fold_change, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, size = 1) +
  geom_point(position = position_jitter(0.15), shape = 21, size = 2, alpha = 0.6) +
  stat_compare_means(method = "wilcox.test", label = "p.format", comparisons = ln_comparisons, size = 4) +
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

# END ####