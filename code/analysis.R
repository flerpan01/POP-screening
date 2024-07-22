# --- pipeline --------------------------------------------------------------- #

# import curated excel file
# log transform the measures
# calculate the mean for each treatment (triplicates) and the fold change
# run mixed linear model, run post-hoc test, save output to excel
# exponent the measures (reverse from log transformation)
# plot

# 230612 make y-axis scales independent for each of the parameters
# (E2, T, Via). Use geom_blank() and "fool" ggplot with extra values.
# Assess the min and max for each compound before plotting. 
# line 223-247 + function make_plot, quick'n'dirty solutions

# --- metadata --------------------------------------------------------------- #

# cmpd = compounds, 24 in total
# treat = treatment, 0 (control) - 3e-05 (max), unit concentration
# n = batch, 1-3

# E2 = estradiol, parameter 1
# T = testosteron, parameter 2
# Via = viability, MTT, parameter 3

# --- packages --------------------------------------------------------------- #

library(openxlsx)
library(tidyverse)
library(lme4)
library(car)
library(cowplot)
library(scales)

library(ggpubr)   # signficance bars

library(multcomp) # post-hoc test
library(broom)    # post-hoc test

# = functions ================================================================ #

make_df <- function(res, treatments, parameter) {
  res <- summary(res)
  
  rows <- 2:length(treatments)
  
  fold_change <- exp(res$coefficients[rows, 1])
  pval <- 2*(1 - pnorm(abs(res$coefficients[rows, 3])))
  upper_ci <- exp(res$coefficients[rows, 1] + 1.96*res$coefficients[rows, 2])
  lower_ci <- exp(res$coefficients[rows, 1] - 1.96*res$coefficients[rows, 2])
  
  data <- data.frame(treatment = as.numeric(as.character(treatments[rows])),
                     fold_change,
                     p_value = pval,
                     upper_ci,
                     lower_ci,
                     parameter = parameter)
  return(data)
}
lm_output <- function(res, parameter) {
  aov <- Anova(res, type = "III")
  out <- data.frame(variable = rownames(aov), aov, parameter = parameter)
  rownames(out) <- NULL
  
  #summary(res)$coefficients
  
  out
}
make_data <- function(d, compound) {
  # subset a compound
  print(compound)
  dat <- subset(d, cmpd %in% compound)
  treatments <- unique(dat$treat)
  
  dat$treat <- factor(dat$treat)
  dat$n <- factor(dat$n)
  
  lm <- list()
  
  lm$Estradiol <- lmer(E2 ~ treat + n + (1 | treat : n), data = dat, REML = TRUE)
  lm$Testosterone <- lmer(T ~ treat + n + (1 | treat : n), data = dat, REML = TRUE)
  lm$Viability <- lmer(Via ~ treat + n + (1 | treat : n), data = dat, REML = TRUE)
  
  # save post-hoc output
  l2 <- list()
  for (parameter in names(lm)) {
    l2[[parameter]] <- data.frame(tidy(glht(lm[[parameter]]))[, -2], parameter = parameter)
  }
  l2 <- Reduce(function(x,y) rbind(x,y), l2)
  l2$compound <- compound
  
  # save anova output
  l1 <- list()
  for (parameter in names(lm)) {
    l1[[parameter]] <- lm_output(lm[[parameter]], parameter)
  }
  output <- Reduce(function(x,y) rbind(x,y), l1)
  output$compound <- compound
  
  # save fold change + CI output
  l1 <- list()
  for (i in names(lm)) {
    l1[[i]] <- make_df(lm[[i]], treatments, i)
  }
  
  data <- Reduce(function(x,y) rbind(x,y), l1)
  rownames(data) <- NULL
  data$compound <- compound
  
  out <- list()
  
  out$data <- data
  out$output <- output
  out$post_hoc <- l2
  
  return(out)
}
make_plot <- function(data, compound, combined = FALSE, facet = "horizontal", pltsize = 14, pointsize = 4, ylim = NULL, sign_ind_size = 2) {
  if (!combined) {


    # table for plotting significance indicators: *, **, ***, ns = nothing
    statdata <- data.frame(
      parameter = data$parameter,
      y.position = data$upper_ci * 1.025, # needed to see sign-indicator
      group1 = data$treatment,
      group2 = data$treatment,
      p_signif = sapply(data$p_value, function(x){
        if (x < 0.001) return("***")
        if (x < 0.01) return("**")
        if (x <= 0.05) return("*")
        if (x > 0.05) return(NA)
      })
    )

    p <- ggplot(data, aes(treatment, fold_change)) +
      geom_hline(yintercept = 1, linetype = 2, linewidth = 0.5, 
        color = "darkorange3") +
      geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +
      geom_line(aes(group = 1), alpha = 0.7, linewidth = 0.4) +
      geom_point(aes(fill = parameter), shape = 21, size = pointsize) +

      {if (facet == "horizontal") facet_grid(parameter ~ .)} +
      {if (facet == "free") facet_grid(parameter ~ ., scales = "free")} +
      {if (facet == "vertical") facet_wrap(. ~ parameter, scales = "free")} +

      scale_x_log10() +
      scale_y_continuous(labels = comma) +
      theme_pubr(pltsize) + 
      theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        strip.text.x = element_text(size = 10)) +

      # set individual y-axis for each compound
      # use geom_blank with max / min values for each compound
      {if (!is.null(ylim)) geom_blank(data = ylim, aes(treatment, fold_change))} +

      labs(x = "Concentrations [M]", y = "Fold of control", title = compound) +

      stat_pvalue_manual(statdata, label = "p_signif", size = sign_ind_size, 
        tip.length = 0)

  }
  
  if (combined) {
    pos <- 0.2
    p <- ggplot(
      data, 
      aes(
        x = treatment, 
        y = fold_change, 
        ymin = lower_ci, 
        ymax = upper_ci, 
        group = parameter, 
        fill = parameter)
    ) +
      
      geom_line(aes(color = parameter), 
                position = position_dodge2(width = pos), linetype = 2) +
      
      geom_errorbar(aes(color = parameter), 
                    width = pos, 
                    position = position_dodge2(width = pos),
                    linewidth = 0.5) +
        
      geom_point(shape = 21, size = 4, position = position_dodge2(width = pos)) +
      scale_y_continuous(breaks = breaks_width(0.5)) +
      scale_x_log10() +
      theme_linedraw() +
      labs(x = "Concentration", y = "Fold Change", title = compound, fill = "ASDASD") +
      guides(color = "none") + 

      theme(legend.position = "bottom",
            plot.title = element_text(hjust = 0.5)) 
  }
  
  return(p)
} 

# = code ===================================================================== #

d <- read.xlsx("doc/steroidogenesis_dataset.xlsx")

# make measures numeric and log transform
cols <- names(d) %in% c("cmpd", "treat", "n")
d[, !cols] <- sapply(d[, !cols], function(x) log2(as.numeric(x) + 0.5))

compounds <- unique(d$cmpd)

# vertical
tmp <- list()
anova_data <- list()
lm_data <- list()
post_hoc_data <- list()

for (compound in compounds) {
  data <- make_data(d, compound)
  tmp[[compound]] <- data$data
  anova_data[[compound]] <- data$data
  lm_data[[compound]] <- data$output
  post_hoc_data[[compound]] <- data$post_hoc
  
  #plt_list[[compound]] <- make_plot(data$data, compound, combined = FALSE, facet = "vertical", pltsize = 8, pointsize = 2) # FIX HERE
}

anova_data <- Reduce(function(x,y) rbind(x,y), anova_data)
lm_data <- Reduce(function(x,y) rbind(x,y), lm_data)
post_hoc_data <- Reduce(function(x,y) rbind(x,y), post_hoc_data)

# --- stats to excel --------------------------------------------------------- #

write.xlsx(
  list(
    "anova_summary" = anova_data, 
    "lm_summary" = lm_data,
    "post_hoc_summary" = post_hoc_data
  ), 
  file = "data/fit_linear_mixed_model_output.xlsx"
)

# --- save to pdf ------------------------------------------------------------ #

# anova_data contains all min and max fold changes
ylim <- anova_data %>%
  group_by(parameter) %>%
  summarise(
    min = round(min(lower_ci), 1),
    max = round(max(upper_ci) * 1.03, 1) # 3% increase y-axis so sign-ind can fit
  )

# fix dataframe to contain fold_change and treatment (y , x)
ylim <- reshape2::melt(ylim, id.vars = "parameter")
colnames(ylim)[3] <- "fold_change"
ylim$treatment <- anova_data$treatment[1]
ylim$lower_ci <- 1
ylim$upper_ci <- 1

# base p-values on compound-level
dat <- lm_data[lm_data$variable == "treat" & lm_data$Pr..Chisq. <= 0.05, ]

rows <- unique(post_hoc_data$contrast)
rows <- rows[grepl("treat", rows)]
post_hoc_data <- subset(post_hoc_data, contrast %in% rows)

# keep track of which are significant
data <- lapply(tmp, function(x) cbind(x, sign = FALSE))
for (i in seq_len(nrow(dat))){
  comp <- dat$compound[i]
  param <- dat$parameter[i]

  out <- data[[comp]]
  
  contrast <- subset(post_hoc_data, parameter %in% param & compound %in% comp & adj.p.value <= 0.05)$contrast
  contrast <- sub("treat", "", contrast)

  rows <- out$treatment %in% contrast & out$parameter %in% param & out$compound %in% comp

  out$sign[rows] <- TRUE

  data[[comp]] <- out
}

# Set p-value of non-significant post hoc correlations = 1, quick-fix
for (i in seq_len(length(data))){
  sign <- data[[i]]$sign
  data[[i]]$p_value[!sign] <- 1
}

# plot
n <- seq_along(unique(d$cmpd))
per_page <- 4

plt_list <- list() # FIX HERE
for (compound in compounds) {
  plt_list[[compound]] <- make_plot(
    data = data[[compound]], 
    compound, 
    combined = FALSE, 
    facet = "vertical", 
    pltsize = 9, 
    pointsize = 2.5,
    ylim = ylim,
    sign_ind_size = 3
  )
} # FIX THIS AND ABOVE

li <- split(n, ceiling(n / per_page))

pdf("img/plots.pdf", width = 7, height = 9)
for (i in li) print(plot_grid(plotlist = plt_list[i], ncol = 1))
dev.off()