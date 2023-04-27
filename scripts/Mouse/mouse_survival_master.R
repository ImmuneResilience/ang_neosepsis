########################
# Vascular Resilience
# Mouse: Survival Analysis
#######################

# In this script, we perform kaplan-meier analysis generate survival curves:
# CS vs AA
# CS vs LNAME
# CS vs A1+LNAME
# CS vs A1 vs A2
# CS vs AntA2
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
survdata <- read.csv('../../Rdata/Mouse/survival_data.csv')
source('helper_functions/monitoring_helper_functions.R')

# Figure Theme
survival_theme <- theme(
  legend.text = element_text(size = 12)
  ) 

survdata$group %>% unique()
survdata$group <- factor(
  survdata$group, 
  c('CS','Ang1','Ang2','CS + AA',
    'L-NAME','Ang1+L-NAME','Ang1+L-Arg',
    'anti-Ang2','L-Arg','Ang1+Ang2')
)

# Set colors 

dark2 <- brewer.pal(n = 8, name = 'Dark2')[c(2,1,3,4)]
cols <- c("#969696", dark2)

lines <- c(2, 1, 4, 3)

# A1_A2_Saline ------------------------------------------------------------
ang1_ang2 = filter(survdata, experiment_type == "ang1_ang2")
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
antiAng2 = filter(survdata, experiment_type == "anti-Ang2")

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
lname = filter(survdata, experiment_type == "L-NAME")

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
r_lname_ang = filter(survdata, experiment_type %in% c("L-NAME_Ang1_2","L-NAME_Ang1_3"))

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
ang1_larg_treatment = filter(survdata, experiment_type == "L-Arg_Ang1_treatment")

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

# CSAA_CS -----------------------------------------------------------------
csaa_cs <- filter(survdata, experiment_type == "ARA_survival")

custom_survplot(csaa_cs,
                pal = dark2,
                lines = lines,
                custom_theme = survival_theme,
                pval_xy = c(11,6))

ggsave(paste0(output.dir,'mouse.survival.CSAA_CS.pdf'),
       width = 8.5,
       height = 8,
       units = 'cm')

# A1A2_Saline -------------------------------------------------------------
ang1ang2_cs = filter(survdata, experiment_type == "Ang1-Ang2 Survival")

custom_survplot(ang1ang2_cs,
                pal = dark2,
                lines = lines,
                custom_theme = survival_theme,
                pval_xy = c(12,5))

ggsave(paste0(output.dir,'mouse.survival.A1A2_Saline.pdf'),
       width = 8.5,
       height = 8,
       units = 'cm')
