####################
# Vascular resilience
# Liver ROS quantification: mouse
###################

# In this script, we produce box plots for ROS quantification of AA and L-NAME treated mice, and compare groups using the wilcoxon test.

# Setup -------------------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load packages
library(plyr)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)

# load data
ros <- read.csv("../../Rdata/Mouse/ms_ros.csv")

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

# LNAME vs CS vs Control --------------------------------------------------

# for this comparison, we use the 2nd batch, and select lname_cs, cs, and control variable.

# subset data
ln <- filter(ros, treatment %in% c("lname_cs", "cs", "control"), batch == 2)

# prepare variables
ln$treatment <- gsub("lname_cs", "lname+cs", ln$treatment)

# prepare comparisons

ln_comparisons <-  list(c("control", "cs"),
                        c("cs", "lname+cs"),
                        c("control", "lname+cs"))

# set factor levels
ln$treatment <- factor(ln$treatment, levels = c("control", "cs", "lname+cs"))

# plot

ggplot(ln, aes(treatment, fold_change, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, size = 1) +
  geom_point(position = position_jitter(0.15), shape = 21, size = 2, alpha = 0.6) +
  stat_compare_means(method = "wilcox.test", label = "p.format", comparisons = ln_comparisons, size = 4) +
  scale_fill_manual(values = cols) +
  ggtitle("LNAME: Liver ROS") +
  theme_classic() + 
  box.theme +
  ylab("Fold Change")

# save plot
ggsave(paste0(output.dir, "mouse_ROS_lname_cs_control.pdf"),
       device = "pdf",
       dpi = 300,
       width = 3, height = 3.5)



# AA vs CS cs Control -----------------------------------------------------

# For this analysis, we use data across all batches, but select only aa, cs, and control treatments

aa <- filter(ros, treatment %in% c("aa_cs", "control", "cs"))
table(aa$treatment)

# set variable names
aa$treatment <- gsub("aa_cs", "aa+cs", aa$treatment)

# set comparisons

aa_comparisons = list(c("aa+cs", "cs"),
                      c("control", "cs"),
                      c("aa+cs", "control"))

# set factor levels
aa$treatment <- factor(aa$treatment, levels = c("control", "cs", "aa+cs"))

# plot
ggplot(aa, aes(treatment, fold_change, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, size = 1) +
  geom_point(position = position_jitter(0.15), shape = 21, size = 2, alpha = 0.6) +
  stat_compare_means(method = "wilcox.test", label = "p.format", comparisons = aa_comparisons, size = 4) +
  scale_fill_manual(values = cols) +
  ggtitle("AA: Liver ROS") +
  theme_classic() + 
  box.theme +
  ylab("Fold Change")

# save plot
ggsave(paste0(output.dir, "mouse_ROS_aa_cs_control.pdf"),
       device = "pdf",
       dpi = 300,
       width = 3, height = 3.5)
# END ####