hours_alive <- function(df){
  # I think easiest way to do this is to melt on HPC columns and pull the value where they are recorded dead
  
  cols_hpc <- colnames(df)[grepl("hpc", colnames(df))]
  cols_status <- colnames(df)[grepl("status", colnames(df))]
  
  dmelt <- reshape2::melt(df, measure.vars = cols_hpc) %>%
    dplyr::rename(hpc = value) %>%
    mutate(visit = substr(variable,0,2), variable = NULL)
    
  dmelt <- reshape2::melt(dmelt, measure.vars = cols_status) %>%
    dplyr::rename(status = value) %>%
    filter(substr(variable,0,2)== visit)
    
  dmelt <- dmelt %>% mutate(
    simple_status = ifelse(is.na(status), NA, 
                           ifelse(status == "EE", "EE", "die"))) 
  
  dmelt <- dmelt %>% filter(is.na(simple_status) == F) %>% 
    dplyr::rename(final_visit = visit, hours_alive = hpc)
  out <- merge(df, dmelt[,c("pup_id", "final_visit", "hours_alive")], by = "pup_id")

  return(out)
}

binary_outcome <- function(x){
  x <- x %>% mutate(death = ifelse(is.na(outcome), NA, ifelse(outcome == "EE", 0, 1)))
}

custom_survplot <- function(Data, pal, lines, max_hour = 72, break_by = 12, pval_xy = c(10,15), pval_size = 4, custom_theme = theme()){
  Data$group <- factor(Data$group)
  # Perform the log-rank test and get p-value
  test <- survdiff(Surv(sac_hpc, outcome) ~ group, data  = Data, rho = 0)
  pval <- broom::glance(test)$p.value 
  # p-value annotation
  anot <- paste0("p = ", round(pval, digits = 5))
  plot <- ggsurvplot(survfit(Surv(sac_hpc, outcome) ~ group, data = Data), 
                     data = Data, 
                     palette = pal,
                     linetype = lines,
                     censor = FALSE, fun = 'pct', 
                     xlab = "Hours post challenge", ylab = "Percent survival",
                     size = 1, 
                     break.x.by = break_by,
                     xlim = c(0,max_hour),
                     risk.table = F,
                     pval = F,
                     legend = "top",
                     legend.labs = names(table(Data$group)),
                     legend.title = element_blank(),
                     tables.theme = theme(title = element_blank())
  )
  
  plot$plot <- plot$plot + annotate("text", size = pval_size, x = pval_xy[1], y = pval_xy[2], label = anot)
  plot$plot <- plot$plot + custom_theme
  plot$plot <- plot$plot + labs(color  = "", linetype = "")
  return(plot)
}
