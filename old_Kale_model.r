rm(list=ls())
library(dplyr)
library(RTMB)
library(tidyr)
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
  group_by(t_k,r_k,t_l, r_l, tag, t_yr, f_yr, f_tk, f_tl) %>%
  summarise(n = n()) %>%
  filter((t_l<=r_l) | is.na(r_l))

formulas <- list(
  phi     = ~ 1 ,#+ (1|f_yr),   # survival
  p       = ~ 1 ,#+ (1|f_yr),   # detection
  v       = ~ 1, #Tagging rate
  w       = ~ 1 + (1|f_tk), #entry
  Nsuper  = ~ 1 ,#+ (1|f_yr),   # Super pop
  t_var   = ~ 1 #transition matrix
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

u <- rep(0,max(d$t_k))
u_i <- aggregate(list(u = d$n),by = list(s = d$t_k),sum)
data$u <- u[u_i$s] <- u_i$u

sp_count <- d %>%
  group_by(t_l) %>%
  summarise(n = sum(n))
data$sp_pr <- as.vector(sp_count$n/sum(sp_count$n))


# make_param_list('~1', d, )
params <- make_all_params(formulas
                          , d
                          , s = max(na.omit(c(d$t_k,d$r_k)))
                          , state_col = data$state)


# make_param_list('~1', d, )
struct <- extract_all_structures(formulas
                                 , d
                                 , s = max(na.omit(c(d$t_k,d$r_k)))
                                 , state = data$state)

data$struct <- struct
data$sp_count <- sp_count

Nsuper_g <- data$struct$groups$Nsuper[[1]]
if(length(Nsuper_g)>1){
  data$uTot_g <- as.numeric( with(data, tapply(n, Nsuper_g, sum)) )
}else{
  data$uTot_g <- sum(data$uTot)
}


n_state <- if(!is.null(data$state)) length(unique(data[[data$state]])) else 1L

map <- list()
random <- c( "logit_v")
if(length(data$struct$Z_list$Nsuper)>0){
  random <- c(random, "u_Nsuper")
}
if(length(data$struct$Z_list$phi)>0){
  random <- c(random, "u_phi")
}
if(length(data$struct$Z_list$p)>0){
  random <- c(random, "u_p")
}
if(length(data$struct$Z_list$w)>0){
  random <- c(random, "u_w")
  map$log_sd_w <- as.factor(NA)
}
if(length(params$t_var)>1){
  random <- c(random, "t_var")
}


environment(model) <- .GlobalEnv
map <- list(M_sigma = as.factor(NA),
            log_sd_w = as.factor(NA)
            # logit_v = as.factor(rep(NA,length(params$logit_v)))
            )
obj <- RTMB::MakeADFun(func = model,
                       parameters = params,
                       data = data,
                       random = random,
                       map = map,
                       silent = FALSE)

# Optimize Model
opt <- nlminb(obj$par, obj$fn, obj$gr)
# sdr <- sdreport(obj)
rep <- obj$report()


