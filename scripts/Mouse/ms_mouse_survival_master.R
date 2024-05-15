########################
# Vascular Resilience
# Mouse: Survival Analysis
#######################

# In this script, we perform kaplan-meier analysis generate survival curves:
# CS vs AA
# CS vs LNAME
# CS vs A1+LNAME
# CS vs A1 vs A2
# CS vs AntiA2
# CS vs A1+A2

# Setup  --------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#set output
output.dir <- "../../Figures/mouse/survival_curves/"

# Load packages
library(plyr)
library(tidyverse)
library(survival)
library(survminer)
library(RColorBrewer)

# load data
survdata <- read.csv('../../Rdata/Mouse/ms_survival_data.csv')
source('helper_functions/monitoring_helper_functions.R')

# Figure Theme
survival_theme <- theme(
  legend.text = element_text(size = 12)
) 

survdata$group <- ifelse(survdata$group == "CS + AA", "AA", as.character(survdata$group))
unique(survdata$group)

survdata$group %>% unique()
survdata$group <- factor(
  survdata$group, 
  c('CS','Ang1','Ang2','AA',
    'L-NAME','Ang1+L-NAME','Ang1+L-Arg',
    'anti-Ang2','L-Arg','Ang1+Ang2')
)

# Set colors 

dark2 <- brewer.pal(n = 8, name = 'Dark2')[c(2,1,3,4)]
cols <- c("#969696", dark2)

lines <- c(2, 1, 4, 3)

# A1_A2_Saline ------------------------------------------------------------
ang1_ang2 <- survdata %>%
  filter(experiment_type == "ang1_ang2") %>%
  droplevels()

table(ang1_ang2$group)

p <- custom_survplot(ang1_ang2,
                     pal = dark2,
                     lines = lines,
                     custom_theme = survival_theme,
                     pval_xy = c(12, 5))

p

ggsave(paste0(output.dir,'mouse.survival.A1_A2_Saline.pdf'),
       width = 8.5,
       height = 8,
       units = 'cm')

# antiA2_Saline -----------------------------------------------------------
antiAng2 <- survdata %>%
  filter(experiment_type == "anti-Ang2") %>%
  droplevels()

table(antiAng2$group)

custom_survplot(antiAng2,
                pal = dark2,
                lines = lines,
                custom_theme = survival_theme,
                pval_xy = c(11,5))


ggsave(paste0(output.dir,'mouse.survival.antiA2_Saline.pdf'),
       width = 8.5,
       height = 8,
       units = 'cm')

# LNAME_Saline ------------------------------------------------------------
lname <- survdata %>%
  filter(experiment_type == "L-NAME") %>%
  droplevels()

table(lname$group)

custom_survplot(lname,
                pal = dark2,
                lines = lines,
                custom_theme = survival_theme,
                pval_xy = c(55,100))

ggsave(paste0(output.dir,'mouse.survival.LNAME_Saline.pdf'),
       width = 8.5,
       height = 8,
       units = 'cm')

# A1+LNAME_LNAME ----------------------------------------------------------
r_lname_ang <- survdata %>%
  filter(experiment_type == "L-NAME_Ang1") %>%
  droplevels()

table(r_lname_ang$group)

custom_survplot(r_lname_ang,
                pal = dark2[c(2:3)],
                lines = lines,
                custom_theme = survival_theme,
                pval_xy = c(11,5))

ggsave(paste0(output.dir,'mouse.survival.A1+LNAME_LNAME.pdf'),
       width = 8.5,
       height = 8,
       units = 'cm')

# A1LArg_A1_LArg_Saline ---------------------------------------------------
ang1_larg_treatment <- survdata %>%
  filter(experiment_type == "L-Arg_Ang1_treatment") %>%
  droplevels()

table(ang1_larg_treatment$group)

custom_survplot(ang1_larg_treatment,
                pal = dark2,
                lines = lines,
                custom_theme = survival_theme,
                pval_xy = c(6,6),
                pval_size = 5)


ggsave(paste0(output.dir,'mouse.survival.A1LArg_A1_Larg_Saline.pdf'),
       width = 11,
       height = 8,
       units = 'cm')

# AA_CS -----------------------------------------------------------------
csaa_cs <- survdata %>%
  filter(experiment_type == "ARA_survival") %>%
  droplevels()

table(csaa_cs$group)

custom_survplot(csaa_cs,
                pal = dark2,
                lines = lines,
                custom_theme = survival_theme,
                pval_xy = c(11,6))

ggsave(paste0(output.dir,'mouse.survival.AA_CS.pdf'),
       width = 8.5,
       height = 8,
       units = 'cm')

# A1A2_Saline -------------------------------------------------------------
ang1ang2_cs <- survdata %>%
  filter(experiment_type == "Ang1-Ang2 Survival") %>%
  droplevels()

table(ang1ang2_cs$group)

custom_survplot(ang1ang2_cs,
                pal = dark2,
                lines = lines,
                custom_theme = survival_theme,
                pval_xy = c(12,5))

ggsave(paste0(output.dir,'mouse.survival.A1A2_Saline.pdf'),
       width = 8.5,
       height = 8,
       units = 'cm')
# END ####