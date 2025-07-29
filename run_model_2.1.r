rm(list=ls())
library(dplyr)
library(RTMB)
library(tidyr)
library(rlang)


frms <- list(#frm1 = list(
  # phi    = ~ 1,
  # p      = ~ 1,
  # w      = ~ -1 + (1|f_tk),
  # t_var = ~ 1,
  # Nsuper = ~ 1)
  # ,frm2 =   list(
  #   phi    = ~ -1 + f_tk,
  #   p      = ~ -1 + f_tk,
  #   w      = ~ -1 + (1|f_tk),
  #   t_var = ~ 1,
  #   Nsuper = ~ 1)
  frm3 =   list(
    phi    = ~ 1 + (1|f_tk),
    p      = ~ 1 + (1|f_tk),
    w      = ~ -1 + (1|f_tk),
    t_var = ~ 1,
    Nsuper = ~ 1)
  # ,frm4 =   list(
  #   phi    = ~ -1 + f_sex + (1|f_tk),
  #   p      = ~ -1 + f_sex + (1|f_tk),
  #   w      = ~ -1 + (1|f_tk),
  #   t_var = ~ 1,
  #   Nsuper = ~ -1 + f_sex)
  # ,frm5 =   list(
  #   phi    = ~ 1 + f_sex + f_tl + (1|f_tk),
  #   p      = ~ 1 + f_sex + (1|f_tk),
  #   w      = ~ -1 + (1|f_tk),
  #   t_var = ~ 1,
  #   Nsuper = ~ 1 + f_sex + f_tl)
)

dans <- read.csv("data/NF Lewis_Chinook_abundance by year and run.csv")
years <- 2024
Est <- expand.grid(form = names(frms), years = years, val = 0)
rep <- list()

for(i in 1:nrow(Est)){
  d <- read.csv("C:/R_packages/kale/data/NF_Lewis_combined-2013_2024-2025-04-18.csv") %>%
    mutate(t_k = Period_Cap,
           r_k = Period_Recap,
           t_yr = Return_Yr) %>%
    filter(Return_Yr == Est$years[i]) %>%
    # mutate(t_l = if_else(is.na(TagState),NA,1),
    #        r_l = if_else(is.na(RecapState),NA,1)) %>%
    # mutate(t_l = TagState,
    #        r_l = RecapState) %>%
    mutate(t_l = if_else(is.na(TagState),NA,TagState),
           r_l = if_else(is.na(RecapState),NA,TagState)) %>%
    # filter(is.na(r_wk) | r_k>0) %>%
    mutate(tag = ifelse(is.na(Tag1),FALSE,TRUE)) %>%
    mutate(f_yr = factor(t_yr, levels = unique(sort(as.numeric(t_yr))))) %>%
    mutate(f_tk = factor(t_k, levels = unique(sort(as.numeric(t_k))))) %>%
    mutate(f_tl = factor(t_l, levels = unique(sort(as.numeric(t_l))))) %>%
    mutate(f_rl = factor(r_l, levels = unique(sort(as.numeric(t_l))))) %>%
    mutate(f_FL = factor(FL, levels = unique(sort(as.numeric(FL))))) %>%
    mutate(f_sex = factor(Sex_final, levels = unique(Sex_final))) %>%
    # filter(f_sex != "Jack") %>%
    group_by(t_k,r_k,t_l, r_l, tag, t_yr, f_yr, f_tk, f_rl, f_tl, f_sex) %>%
    summarise(n = n()) %>%
    filter((t_l<=r_l) | is.na(r_l))



  formulas <- frms[[Est$form[i]]]

  input <- list(state = NULL,
                time = "f_tk")

  source("run_tmb.r")

  Est$val[i] <- sum(out$rep$Nsuper)
  rep[[i]] <- out$rep
}
Est$obs <- dans$est[match(Est$years,dans$Year)]
print(Est)
