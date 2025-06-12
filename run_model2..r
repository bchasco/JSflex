rm(list=ls())
library(dplyr)
library(RTMB)
library(tidyr)
d <- read.csv("C:/R_packages/kale/data/NF_Lewis_combined-2013_2024-2025-04-18.csv") %>%
  mutate(t_k = Period_Cap,
         r_k = Period_Recap,
         t_yr = Return_Yr) %>%
  filter(Return_Yr == 2020) %>%
  mutate(t_l = TagState,
         r_l = RecapState) %>%
  # filter(is.na(r_wk) | r_k>0) %>%
  mutate(tag = ifelse(is.na(Tag1),FALSE,TRUE)) %>%
  mutate(f_yr = factor(t_yr, levels = unique(sort(as.numeric(t_yr))))) %>%
  mutate(f_tk = factor(t_k, levels = unique(sort(as.numeric(t_k))))) %>%
  mutate(f_tl = factor(t_l, levels = unique(sort(as.numeric(t_l))))) %>%
  mutate(f_sex = factor(Sex_final, levels = unique(Sex_final))) %>%
  group_by(t_k,r_k,t_l, r_l, tag, t_yr, f_yr, f_tk, f_tl, f_sex) %>%
  summarise(n = n()) %>%
  filter((t_l<=r_l) | is.na(r_l))

formulas <- list(
  phi     = ~ 1 ,   # survival
  p       = ~ 1 ,   # detection
  v       = ~ 1 , #Tagging rate
  w       = ~ 1 , #entry
  Nsuper  = ~ 1 ,   # Super pop
  t_var   = ~ 1#transition matrix
)


data <- list(
  s = max(d$t_k, na.omit(d$r_k)),                    # Number of periods
  t_k = d$t_k,               # First detection
  r_k = d$r_k,               # Recapture time (NA if none)
  t_l = d$t_l,               # First detection
  r_l = d$r_l,               # Recapture time (NA if none)
  tag = as.integer(d$tag),   # 1 = tagged, 0 = untagged
  n = d$n,                   # Number of fish in each group
  uTot = sum(d$n),           # Total number of detections
  sim_mode = FALSE,
  # state = NULL
  state = "t_l"
)

source("model_source.r")
source("oldPlot.r")
print(rep$out$Nsuper_vec)
print(sum(rep$out$Nsuper_vec))
