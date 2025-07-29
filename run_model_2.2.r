rm(list=ls())
library(dplyr)
library(RTMB)
library(tidyr)
library(rlang)
d <- read.csv("C:/R_packages/kale/data/NF_Lewis_combined-2013_2024-2025-04-18.csv") %>%
  mutate(t_k = Period_Cap,
         r_k = Period_Recap,
         t_yr = Return_Yr) %>%
  filter(Return_Yr == 2024) %>%
  mutate(t_l = TagState,
         r_l = RecapState) %>%
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

# Survival ~ (1|sex) + (1|reach) + (1|time)
# Capture ~ (1|sex) + (1|reach) + (1|time)
# Entry ~ (1|reach/sex) + (1|time/sex)
# Tagging rate ~ I don’t care
#
#
# Some further enhancements:
#   Short Term:
#   1)	I’d like to try to run the models at a daily time step since “periods” may not actually mean the same thing across space (if sampling dates differ among reaches)
# Medium Term:
#   2)	Ideally rather than simple random effects, time and space would be smoothers (e.g., spline or random walk for time, OU for space)
# 3)	We could also look at ways to fit spatiotemporal random effects in datasets that will support it
# Longer term?:
#   4)	Incorporating individual continuous covariates could be cool…e.g., length
# 5)	I’d like to eventually explore GPSing tag deployment and recovery locations so that we have the flexibility to “create” our reaches as part of the analysis…and then can explore the sensitivity of the estimates to reaches of arbitrary size


##-----------------------------------------------
## 1.1  design matrices and flists  -------------
##-----------------------------------------------
formulas <- list(
  phi    = ~ 1 + f_sex,
  p      = ~ 1 + (1|f_tl),
  w      = ~ 1,
  v = ~ 1 , #deprecated for now
  t_var = ~ 1,
  Nsuper = ~ 1
)

input <- list(
  state = "f_tl",
  period = "f_tk"
)


design <- make_design_list_2(formulas, input, as.data.frame(d))

data <- make_RTMB_data_list(design, d, state = design$state, period = "f_tk")
params <- make_param_vectors(design, d)
random <- make_random(params)

G <- max(lk$phi$lookup$g)
S <- if (!is.null(design$state)) nlevels(d$f_tl) + 2 else 3L
T <- max(lk$phi$lookup$time)

S2 <- array(0, c(G, S, S, T))

map <- make_map(params)
#
#
source("model_2.2.r")
obj <- RTMB::MakeADFun(model,
                       data = data,
                       random = random,
                       map = map,
                       parameters = params)
# # Optimize Model
# opt <- nlminb(obj$par, obj$fn, obj$gr)
rep <- obj$report()
# print(dim(rep$S2))
# print(length(params$u_t_var))
# # sdr <- sdreport(obj)
# #
# # print(dim(rep$xz$phi$Z))
# # print(dim(rep$xz$phi$X))
# # print(dim(rep$phi_combo_df))
# #
plot_Sankey(data)
# #
