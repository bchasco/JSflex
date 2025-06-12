# d_raw <- read.csv("C:/R_packages/kale/data/NF_Lewis_combined-2013_2024-2025-04-18.csv") %>%
#   mutate(t_k = Period_Cap,
#          r_k = Period_Recap,
#          t_yr = Return_Yr) %>%
#   filter(Return_Yr > 2018) %>%
#   mutate(t_l = TagState,
#          r_l = RecapState) %>%
#   # filter(is.na(r_wk) | r_k>0) %>%
#   mutate(tag = ifelse(is.na(Tag1),FALSE,TRUE)) %>%
#   mutate(f_yr = factor(t_yr, levels = unique(sort(as.numeric(t_yr))))) %>%
#   mutate(f_tk = factor(t_k, levels = unique(sort(as.numeric(t_k))))) %>%
#   group_by(t_k,r_k,t_l, r_l, tag, t_yr, f_yr, f_tk) %>%
#   summarise(n = n()) %>%
#   filter((t_l<=r_l) | is.na(r_l))
#
# formulas <- list(
#   phi = ~ 1 + (1|f_yr),
#   p = ~ 1 + (1|f_yr)  ,
#   psi = ~ 1,
#   v = ~ 1,
#   w = ~ 1 ,
#   Nsuper = ~ 1 + (1|f_yr)
# )
#
# fit <- fit_jsflex(d_raw,
#                   formulas = formulas)
#

fit()
