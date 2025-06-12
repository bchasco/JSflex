
library(tidyr)
library(dplyr)
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
  group_by(t_k,r_k,t_l, r_l, tag, t_yr, f_yr) %>%
  summarise(n = n()) %>%
  filter((t_l<=r_l) | is.na(r_l))

u <- rep(0,max(d$t_k))
u_i <- aggregate(list(u = d$n),by = list(s = d$t_k),sum)
u[u_i$s] <- u_i$u

sp_count <- d %>%
  group_by(t_l) %>%
  summarise(n = sum(n))


data <- list(
  s = max(d$t_k, na.omit(d$r_k)),                    # Number of periods
  t_k = d$t_k,               # First detection
  r_k = d$r_k,               # Recapture time (NA if none)
  t_l = d$t_l,               # First detection
  r_l = d$r_l,               # Recapture time (NA if none)
  tag = as.integer(d$tag),   # 1 = tagged, 0 = untagged
  n = d$n,                   # Number of fish in each group
  u = u,          # First detections per period (for multinomial)
  uTot = sum(d$n),           # Total number of detections
  mod = 1, #3 forward matrix (1D), 4 forward matrix (multi D)
  sp_pr = as.vector(sp_count$n/sum(sp_count$n)),
  sim_mode = FALSE
)


formulas <- list(
  phi = ~ 1,
  p = ~ 1,
  psi = ~ 1,
  v = ~ 1,
  w = ~1,
  Nsuper = ~ 1
)

# make_param_list('~1', d, )
params <- make_all_params(formulas
                          , d
                          , s = max(na.omit(c(d$t_k,d$r_k)))
                          , state_col = NULL)


# make_param_list('~1', d, )
struct <- extract_all_structures(formulas
                          , d
                          , s = max(na.omit(c(d$t_k,d$r_k)))
                          , state = NULL)

data$struct <- struct

phi_g <- data$struct$groups$phi$group_phi_re1
if(length(phi_g)>1){
  data$uTot_g <- as.numeric( with(data, tapply(n, phi_g, sum)) )
}else{
  data$uTot_g <- sum(data$uTot)
}
