sex_loc_ttt_model <- function(parms) {

  ## ------------------------------------------------------------------
  ## 0.  unpack data & parameters  ------------------------------------
  ## ------------------------------------------------------------------
  RTMB::getAll(data, parms)

  nll <- 0

  n_loc    <- 5

  ## ------------------------------------------------------------------
  ## 1. transition array  S (G × … × …) ---------------------------------
  T_mat <- array(0, c(G, 11, 11, s))
  for(g in 1:G){
    for (t in seq_len(s)) {
      icnt <- 1
      for(loc in 1:4){
        for(loc2 in (loc+1):5){
          T_mat[g,loc,loc2,t] <- exp(t_vars[1,icnt])
          icnt <- icnt + 1
        }
        T_mat[g,loc,loc,t] <- -sum(T_mat[g,loc,,t])
      }
    }
  }
  RTMB::REPORT(T_mat)
  RTMB::REPORT(t_vars)

  nll_u_t_vars <- 0
  nll_u_t_vars <- nll_u_t_vars -sum(dnorm(t_vars,0,exp(log_sd_t_vars),TRUE))


  S <- array(0, c(G, 11, 11, s))
  g_u_phi <- c(1,2,3,4,5)
  for(g in 1:G){
    for (t in seq_len(s)) {
      for(loc in 1:5){
        # eta_phi <- phi[g,loc]
        eta_phi <- phi[1]
        if(t<s){
          eta_phi <- eta_phi + u_phi_gls[g,g_u_phi[loc],t] + u_phi_gl[g,g_u_phi[loc]] + u_phi_g[g]
        }
        S[g,loc,loc,t] <- -exp(eta_phi)
        S[g,loc,,t] <- S[g,loc,,t] + T_mat[g,loc,,t]
        S[g,loc,11,t] <- -sum(S[g,loc,,t])
      }
      S[g,,,t] <- as.matrix(Matrix::expm(S[g,,,t]))
    }
  }
  RTMB::REPORT(S)
  RTMB::REPORT(u_phi_g)
  RTMB::REPORT(u_phi_gl)
  RTMB::REPORT(u_phi_gls)
  nll_u_phi <- 0
  # nll_u_phi <- -sum(dnorm(u_phi_gls,0,exp(log_sd_phi),TRUE))
  nll_u_phi <- nll_u_phi-sum(dnorm(u_phi_gl,0,exp(log_sd_phi),TRUE))
  # nll_u_phi <- nll_u_phi-sum(dnorm(u_phi_g,0,exp(log_sd_phi),TRUE))

  D <- array(0, dim = c(G, 11, 11, s))
  g_u_p <- c(1,2,3,4,5)
  for(g in 1:G){
    for (t in 1:s) {
      for(loc in 1:5){
        mu_p <- p[g,loc]
        # mu_p <- p[1]
        D[g,loc,loc,t] <- -exp(mu_p + u_p_gls[g,g_u_p[loc],t] + u_p_gl[g,loc] + u_p_g[g])
        D[g,loc,loc+5,t] <- -D[g,loc,loc,t]
      }
      D[g,,,t] <- as.matrix(Matrix::expm(D[g,,,t]))
    }
  }
  RTMB::REPORT(D)
  RTMB::REPORT(u_p_gls)
  RTMB::REPORT(u_p_gl)
  RTMB::REPORT(u_p_g)
  nll_u_p <- 0
  # nll_u_p <- -sum(dnorm(u_p_gls,0,exp(log_sd_p),TRUE))
  # nll_u_p <- nll_u_p -sum(dnorm(u_p_gl,0,exp(log_sd_p),TRUE))
  # nll_u_p <- nll_u_p -sum(dnorm(u_p_g,0,exp(log_sd_p),TRUE))


  Nsuper_vec <- exp(Nsuper)
  RTMB::REPORT(Nsuper_vec)

  delta_mat <- array(0, c(G, 5, s))
  # Get u_uniform and eta_w with same length as combo_df$row
  for(g in 1:G){
    for(loc in 1:5){
      u_uniform  <- pnorm(u_w[1,1,])
      delta_vals <- qgamma(u_uniform, shape = exp(log_sd_w), scale = 1)
      delta_mat[g,loc, ] <- delta_vals
      delta_mat[g,loc,] <- delta_mat[g,loc,]/sum((delta_mat[g,loc,]))

    }
  }
  RTMB::REPORT(u_w)
  RTMB::REPORT(delta_mat)
  pent_mat <- delta_mat

  RTMB::REPORT(pent_mat)
  RTMB::ADREPORT(pent_mat)
  #
  # # ## add prior
  nll_u_w <- 0
  nll_u_w <-  - sum(dnorm(u_w, 0, exp(log_sd_w), log = TRUE))
  RTMB::REPORT(nll_u_w)
  #
  ## ------------------------------------------------------------------
  ## 6.  CJS likelihood (per-fish loop)  ------------------------------
  ## ------------------------------------------------------------------
  CJS_nll <- numeric(length(t_k))
  for (i in seq_along(t_k)) {

    entry      <- t_k[i]
    recapture  <- r_k[i]
    is_tagged  <- tag[i]
    count      <- n[i]


    ## seed α-vector
    alpha      <- numeric(11)
    idx_live   <- if (n_loc > 1) t_l[i] else 1L
    recap_live   <- if (n_loc > 1) f_rl[i] else 1L
    # print(recap_live)
    alpha[idx_live] <- 1

    ## transition kernel for this fish
    ## transition kernel for this fish
    Tmat <- S[sex[i],,,t_k[i]] %*% D[sex[i],,,t_k[i]]

    if(entry<s){
      for (t in seq(entry , s - 1)){
        alpha <- alpha %*% Tmat
        Tmat <- S[sex[i],,,t + 1] %*% D[sex[i],,,t + 1]
      }
    }else{
      alpha <- alpha %*% Tmat
    }

    # Set default in case of no recapture
    prob_rec <- 0

    # If recaptured, set detection index and extract from alpha
    if(is_tagged){
      if(!is.na(recapture)){
        prob_rec <- alpha[recap_live + 5]
      }else{
        prob_rec <- sum(alpha[(n_loc + 1):(2*n_loc)])
      }
      # print(paste(i,idx_live,recap_live))
      # print(prob_rec)
    }

    if (entry == s) prob_rec <- 0

    # print(paste(i, is_tagged, recapture))
    ## ------------------- likelihood contributions -------------------
    if (is_tagged) {
      CJS_nll[i] <- - dbinom( ifelse(is.na(recapture), 0, count),
                              count, prob_rec, log = TRUE)
    }
  }
  RTMB::REPORT(CJS_nll)
  # print(CJS_nll)
  ## ------------------------------------------------------------------
  ## 4.  calc_psi (unchanged – S,D,pent now already G-sized) ---------
  ## ------------------------------------------------------------------
  # print(G)
  calc_psi <- function(Smat, Dmat, pent_mat, sp_pr, s, n_loc) {
    psiPtot_g <- matrix(0,G,5); detect_mat <- array(0, dim = c(G, 5, s))

    sp_pr <- 1
    init_t <- 1

    for(g in 1:data$G){
      #Forward algorithm
      for (j in 1:5){
        alpha <- numeric(11);
        alpha[1:5] <- pent_mat[g,,init_t] * sp_pr #First entrants
        for (t in init_t:s) {
          ## entrants during week t
          alpha[j] <- alpha[j] +
            pent_mat[g,j,t] * sp_pr * (Smat[g,j,j,t] - 1) / log(Smat[g,j,j,t])
          p_vec   <- Dmat[g , 1:5 , 1:5 + 5,t]
          detect_mat[g,j,t] <- sum(alpha[1:n_loc] * p_vec)
          alpha <- alpha %*% Smat[g,,,t] %*% Dmat[g,,,t]
        }
        psiPtot_g[g,j] <- sum(detect_mat[g,j,])
      }
    }

    list(psiPtot_g=psiPtot_g, detect_mat=detect_mat)
  }
  #
  out_g <- calc_psi(S, D, pent_mat, sp_pr, s, 5)
  psiPtot_g  <- out_g$psiPtot_g
  RTMB::REPORT(psiPtot_g)

  detect_mat <- out_g$detect_mat
  RTMB::REPORT(detect_mat)

  lambda_k <- Nsuper_vec * psiPtot_g           # length K
  RTMB::REPORT(lambda_k)
  RTMB::REPORT(Nsuper_vec)

  ## ---------------------------------------------------------------
  ## 2. group-wise total detections uTot_g -------------------------
  # uTot_g <- data.frame(n = n, group = sex) %>%
  #   group_by(group) %>%
  #   summarise(uTot = sum(n), .groups = "drop")
  uTot_g <- array(0, dim = c(G,5))
  for (i in seq_along(n)) {
    g <- sex[i]
    j <- t_l[i]
    uTot_g[g, j] <- uTot_g[g, j] + n[i]
  }
  RTMB::REPORT(uTot_g)
  ## make sure it has length K:
  uTot_g_nll <- -sum(dpois(uTot_g, lambda_k, log = TRUE))
  RTMB::REPORT(uTot_g_nll)
  RTMB::ADREPORT(lambda_k)

  print("test")
  ## ---------------------------------------------------------------
  ## 3. week × group counts  u_k_t  --------------------------------
  ## Build k_idx via Nsuper_for_joint
  ## Build observation matrix by joint group g × period t
  u_g_t <- array(0, dim = c(G,5,s))
  for (i in seq_along(n)) {
    t <- t_k[i]
    g <- sex[i]
    j <- t_l[i]
    u_g_t[g, j, t] <- u_g_t[g, j, t] + n[i]
  }
  RTMB::REPORT(u_g_t)

  ## Multinomial likelihood over each joint group g
  uTot_g_t_nll <- 0
  for (g in 1:G) {
    for (j in 1:5) {
      prob <- (detect_mat[g,j, ] + 1e-6)
      prob <- prob / sum(prob)
      uTot_g_t_nll <- uTot_g_t_nll - dmultinom(
        x = u_g_t[g,j, ],
        prob = prob,
        log = TRUE
      )
    }
  }
  RTMB::REPORT(uTot_g_t_nll)
  RTMB::ADREPORT(prob)

  nll <- sum(CJS_nll)
  nll <- nll + nll_u_w + nll_u_phi + nll_u_p + nll_u_t_vars
  nll <- nll + uTot_g_nll + uTot_g_t_nll
  #
  Ntotal <- sum(Nsuper_vec)
  RTMB::ADREPORT(Nsuper_vec)
  RTMB::ADREPORT(Ntotal)

  return(nll)
}

# rm(list=ls())
library(dplyr)
library(RTMB)
library(tidyr)
library(rlang)
#
#
#
# dans <- read.csv("data/NF Lewis_Chinook_abundance by year and run.csv")
years <- 2013:2024
# my_rep <- list()
# my_sd <- list()
#

#ttt
for(yyy in c(7,8)){
  f <- read.csv("C:/R_packages/kale/data/NF_Lewis_combined-2013_2024-2025-04-18.csv")
  # yyy <- 9
    d <- f %>%
      mutate(t_k = Period_Cap,
           r_k = Period_Recap,
           t_yr = Return_Yr) %>%
    filter(Return_Yr == years[yyy]) %>%
    mutate(t_l = if_else(is.na(TagState),NA,TagState),
           r_l = if_else(is.na(RecapState),NA,RecapState)) %>%
    mutate(tag = ifelse(is.na(Tag1),FALSE,TRUE)) %>%
    mutate(f_yr = factor(t_yr, levels = unique(sort(as.numeric(t_yr))))) %>%
    mutate(f_tk = factor(t_k, levels = unique(sort(as.numeric(t_k))))) %>%
    mutate(f_tl = factor(t_l, levels = unique(sort(as.numeric(t_l))))) %>%
    mutate(f_rl = factor(r_l, levels = unique(sort(as.numeric(r_l))))) %>%
    mutate(f_FL = factor(FL, levels = unique(sort(as.numeric(FL))))) %>%
    mutate(f_sex = factor(Sex_final, levels = unique(Sex_final))) %>%
      mutate(f_sex = if_else(f_sex!="Jack","Adult","Jack")) %>%
      mutate(f_sex = factor(f_sex)) %>%
      mutate(p_f_sex = f_sex) %>%
      mutate(p_f_sex = if_else(p_f_sex!="Jack","Adult","Jack")) %>%
      mutate(p_f_sex = factor(p_f_sex)) %>%
      group_by(t_k,r_k,t_l, r_l, tag, t_yr, f_yr, f_tk, f_rl, f_tl, f_sex, p_f_sex) %>%
    summarise(n = n()) %>%
    filter((t_l<=r_l) | is.na(r_l))

  #Number of periods
  s <- max(d$t_k)
  G <- 2

  # Counts per period
  u_vec <- numeric(s)
  u_sum <- aggregate(list(n = d$n), by = list(t = d$t_k), sum)
  u_vec[u_sum$t] <- u_sum$n

  # Spatial proportions
  sp_count <- aggregate(list(n = d$n),
                        by = list(loc = d$f_tl), sum)
  sp_pr <- sp_count$n / sum(sp_count$n)

  uTot_g <- sum(d$n)

  data <- list(
    t_k               = d$t_k,
    r_k               = d$r_k,
    t_l               = d$t_l,
    f_tl              = as.integer(d$f_tl),
    f_rl              = as.integer(d$f_rl),
    r_l               = d$r_l,
    tag               = as.integer(d$tag),
    n                 = d$n,
    s                 = s,
    G = G,
    sex = as.integer(d$f_sex),
    p_sex = as.integer(d$p_f_sex),
    uTot_g            = uTot_g,
    sp_pr             = sp_pr
  )

  params <- list(u_phi_gls = array(0,c(data$G,5,data$s)),
                 u_phi_gl = array(0,c(data$G,5)),
                 u_phi_g = array(0,c(data$G)),
                 log_sd_phi = (-1),
                 u_p_gls = array(0,c(data$G,5,data$s)),
                 u_p_gl = array(0,c(data$G,5)),
                 u_p_g = array(0,c(data$G)),
                 log_sd_p = (-1),
                 u_w = array(0,c(data$G,1,data$s)),
                 log_sd_w = (-1),
                 # phi = array(-0.5,c(data$G,5)),
                 phi = -1,
                 p = array(-0.5,c(data$G,5)),
                 # t_vars = array(0, c(data$G,5*(5-1)/2)),
                 t_vars = array(0, c(1,5*(5-1)/2)),
                 log_sd_t_vars = -1
                 ,Nsuper = array(6.5,c(data$G,5))
  )

  map <- list(u_p_gls = as.factor(array(NA,c(data$G,5,data$s)))
              ,u_p_gl = as.factor(array(NA,c(data$G,5)))
              ,u_p_g = as.factor(array(NA,c(data$G)))
              ,log_sd_p = as.factor(NA)
              ,u_phi_gls = as.factor(array(NA,c(data$G,5,data$s)))
              # ,u_phi_gl = as.factor(array(NA,c(data$G,5)))
              ,u_phi_g = as.factor(array(NA,c(data$G)))
              # ,log_sd_phi = as.factor(NA)
              # ,t_vars = as.factor(rep(NA,length(params$t_vars)))
              # ,log_sd_t_vars = as.factor(NA)
  )



  sex_loc_ttt_obj <- MakeADFun(sex_loc_ttt_model,
                       data = data,
                       random = c("u_phi_gls","u_phi_gl","u_phi_g","u_p_gls","u_p_gl","u_p_g","u_w","t_vars"),
                       map = map,
                       parameters = params)

  opt <- nlminb(sex_loc_ttt_obj$par, sex_loc_ttt_obj$fn, sex_loc_ttt_obj$gr)
  rep <- sex_loc_ttt_obj$report()
  sdr <- sdreport(sex_loc_ttt_obj)
  my_rep[[yyy]] <- rep
  # sdr_ttt <- sdreport(sex_loc_ttt_model)
  my_sd[[yyy]] <- list(est = as.list(sdr, "Estimate", report=TRUE),
                  sd = as.list(sdr, "Std. Error", report=TRUE))
  # sex_loc_ttt_model <- function(parms) {
  #
  #   ## ------------------------------------------------------------------
  #   ## 0.  unpack data & parameters  ------------------------------------
  #   ## ------------------------------------------------------------------
  #   RTMB::getAll(data, parms)
  #
  #   nll <- 0
  #
  #   n_loc    <- 5
  #
  #   ## ------------------------------------------------------------------
  #   ## 1. transition array  S (G × … × …) ---------------------------------
  #   T_mat <- array(0, c(G, 11, 11, s))
  #   for(g in 1:G){
  #     for (t in seq_len(s)) {
  #       icnt <- 1
  #       for(loc in 1:4){
  #         for(loc2 in (loc+1):5){
  #           T_mat[g,loc,loc2,t] <- exp(t_vars[g,icnt])
  #           icnt <- icnt + 1
  #         }
  #         T_mat[g,loc,loc,t] <- -sum(T_mat[g,loc,,t])
  #       }
  #     }
  #   }
  #   RTMB::REPORT(T_mat)
  #   RTMB::REPORT(t_vars)
  #
  #   nll_u_t_vars <- 0
  #   nll_u_t_vars <- nll_u_t_vars -sum(dnorm(t_vars,0,exp(log_sd_t_vars),TRUE))
  #
  #
  #   S <- array(0, c(G, 11, 11, s))
  #   g_u_phi <- c(1,2,3,4,5)
  #   for(g in 1:G){
  #     for (t in seq_len(s)) {
  #       for(loc in 1:5){
  #         # eta_phi <- phi[g,loc]
  #         eta_phi <- phi[1]
  #         if(t<s){
  #           eta_phi <- eta_phi + u_phi_gls[g,g_u_phi[loc],t] + u_phi_gl[g,g_u_phi[loc]] + u_phi_g[g]
  #         }
  #         S[g,loc,loc,t] <- -exp(eta_phi)
  #         S[g,loc,,t] <- S[g,loc,,t] + T_mat[g,loc,,t]
  #         S[g,loc,11,t] <- -sum(S[g,loc,,t])
  #       }
  #       S[g,,,t] <- as.matrix(Matrix::expm(S[g,,,t]))
  #     }
  #   }
  #   RTMB::REPORT(S)
  #   RTMB::REPORT(u_phi_g)
  #   RTMB::REPORT(u_phi_gl)
  #   RTMB::REPORT(u_phi_gls)
  #   nll_u_phi <- 0
  #   # nll_u_phi <- -sum(dnorm(u_phi_gls,0,exp(log_sd_phi),TRUE))
  #   nll_u_phi <- nll_u_phi-sum(dnorm(u_phi_gl,0,exp(log_sd_phi),TRUE))
  #   # nll_u_phi <- nll_u_phi-sum(dnorm(u_phi_g,0,exp(log_sd_phi),TRUE))
  #
  #   D <- array(0, dim = c(G, 11, 11, s))
  #   g_u_p <- c(1,2,3,4,5)
  #   for(g in 1:G){
  #     for (t in 1:s) {
  #       for(loc in 1:5){
  #         mu_p <- p[g,loc]
  #         # mu_p <- p[1]
  #         D[g,loc,loc,t] <- -exp(mu_p + u_p_gls[g,g_u_p[loc],t] + u_p_gl[g,loc] + u_p_g[g])
  #         D[g,loc,loc+5,t] <- -D[g,loc,loc,t]
  #       }
  #       D[g,,,t] <- as.matrix(Matrix::expm(D[g,,,t]))
  #     }
  #   }
  #   RTMB::REPORT(D)
  #   RTMB::REPORT(u_p_gls)
  #   RTMB::REPORT(u_p_gl)
  #   RTMB::REPORT(u_p_g)
  #   nll_u_p <- 0
  #   # nll_u_p <- -sum(dnorm(u_p_gls,0,exp(log_sd_p),TRUE))
  #   # nll_u_p <- nll_u_p -sum(dnorm(u_p_gl,0,exp(log_sd_p),TRUE))
  #   # nll_u_p <- nll_u_p -sum(dnorm(u_p_g,0,exp(log_sd_p),TRUE))
  #
  #
  #   Nsuper_vec <- exp(Nsuper)
  #   RTMB::REPORT(Nsuper_vec)
  #
  #   delta_mat <- array(0, c(G, 5, s))
  #   # Get u_uniform and eta_w with same length as combo_df$row
  #   for(g in 1:G){
  #     for(loc in 1:5){
  #       u_uniform  <- pnorm(u_w[1,1,])
  #       delta_vals <- qgamma(u_uniform, shape = exp(log_sd_w), scale = 1)
  #       delta_mat[g,loc, ] <- delta_vals
  #       delta_mat[g,loc,] <- delta_mat[g,loc,]/sum((delta_mat[g,loc,]))
  #
  #     }
  #   }
  #   RTMB::REPORT(u_w)
  #   RTMB::REPORT(delta_mat)
  #   pent_mat <- delta_mat
  #
  #   RTMB::REPORT(pent_mat)
  #   RTMB::ADREPORT(pent_mat)
  #   #
  #   # # ## add prior
  #   nll_u_w <- 0
  #   nll_u_w <-  - sum(dnorm(u_w, 0, exp(log_sd_w), log = TRUE))
  #   RTMB::REPORT(nll_u_w)
  #   #
  #   ## ------------------------------------------------------------------
  #   ## 6.  CJS likelihood (per-fish loop)  ------------------------------
  #   ## ------------------------------------------------------------------
  #   CJS_nll <- numeric(length(t_k))
  #   for (i in seq_along(t_k)) {
  #
  #     entry      <- t_k[i]
  #     recapture  <- r_k[i]
  #     is_tagged  <- tag[i]
  #     count      <- n[i]
  #
  #
  #     ## seed α-vector
  #     alpha      <- numeric(11)
  #     idx_live   <- if (n_loc > 1) t_l[i] else 1L
  #     recap_live   <- if (n_loc > 1) f_rl[i] else 1L
  #     # print(recap_live)
  #     alpha[idx_live] <- 1
  #
  #     ## transition kernel for this fish
  #     ## transition kernel for this fish
  #     Tmat <- S[sex[i],,,t_k[i]] %*% D[sex[i],,,t_k[i]]
  #
  #     if(entry<s){
  #       for (t in seq(entry , s - 1)){
  #         alpha <- alpha %*% Tmat
  #         Tmat <- S[sex[i],,,t + 1] %*% D[sex[i],,,t + 1]
  #       }
  #     }else{
  #       alpha <- alpha %*% Tmat
  #     }
  #
  #     # Set default in case of no recapture
  #     prob_rec <- 0
  #
  #     # If recaptured, set detection index and extract from alpha
  #     if(is_tagged){
  #       if(!is.na(recapture)){
  #         prob_rec <- alpha[recap_live + 5]
  #       }else{
  #         prob_rec <- sum(alpha[(n_loc + 1):(2*n_loc)])
  #       }
  #       # print(paste(i,idx_live,recap_live))
  #       # print(prob_rec)
  #     }
  #
  #     if (entry == s) prob_rec <- 0
  #
  #     # print(paste(i, is_tagged, recapture))
  #     ## ------------------- likelihood contributions -------------------
  #     if (is_tagged) {
  #       CJS_nll[i] <- - dbinom( ifelse(is.na(recapture), 0, count),
  #                               count, prob_rec, log = TRUE)
  #     }
  #   }
  #   RTMB::REPORT(CJS_nll)
  #   # print(CJS_nll)
  #   ## ------------------------------------------------------------------
  #   ## 4.  calc_psi (unchanged – S,D,pent now already G-sized) ---------
  #   ## ------------------------------------------------------------------
  #   # print(G)
  #   calc_psi <- function(Smat, Dmat, pent_mat, sp_pr, s, n_loc) {
  #     psiPtot_g <- matrix(0,G,5); detect_mat <- array(0, dim = c(G, 5, s))
  #
  #     sp_pr <- 1
  #     init_t <- 1
  #
  #     for(g in 1:data$G){
  #       #Forward algorithm
  #       for (j in 1:5){
  #         alpha <- numeric(11);
  #         alpha[1:5] <- pent_mat[g,,init_t] * sp_pr #First entrants
  #         for (t in init_t:s) {
  #           ## entrants during week t
  #           alpha[j] <- alpha[j] +
  #             pent_mat[g,j,t] * sp_pr * (Smat[g,j,j,t] - 1) / log(Smat[g,j,j,t])
  #           p_vec   <- Dmat[g , 1:5 , 1:5 + 5,t]
  #           detect_mat[g,j,t] <- sum(alpha[1:n_loc] * p_vec)
  #           alpha <- alpha %*% Smat[g,,,t] %*% Dmat[g,,,t]
  #         }
  #         psiPtot_g[g,j] <- sum(detect_mat[g,j,])
  #       }
  #     }
  #
  #     list(psiPtot_g=psiPtot_g, detect_mat=detect_mat)
  #   }
  #   #
  #   out_g <- calc_psi(S, D, pent_mat, sp_pr, s, 5)
  #   psiPtot_g  <- out_g$psiPtot_g
  #   RTMB::REPORT(psiPtot_g)
  #
  #   detect_mat <- out_g$detect_mat
  #   RTMB::REPORT(detect_mat)
  #
  #   lambda_k <- Nsuper_vec * psiPtot_g           # length K
  #   RTMB::REPORT(lambda_k)
  #   RTMB::REPORT(Nsuper_vec)
  #
  #   ## ---------------------------------------------------------------
  #   ## 2. group-wise total detections uTot_g -------------------------
  #   # uTot_g <- data.frame(n = n, group = sex) %>%
  #   #   group_by(group) %>%
  #   #   summarise(uTot = sum(n), .groups = "drop")
  #   uTot_g <- array(0, dim = c(G,5))
  #   for (i in seq_along(n)) {
  #     g <- sex[i]
  #     j <- t_l[i]
  #     uTot_g[g, j] <- uTot_g[g, j] + n[i]
  #   }
  #   RTMB::REPORT(uTot_g)
  #   ## make sure it has length K:
  #   uTot_g_nll <- -sum(dpois(uTot_g, lambda_k, log = TRUE))
  #   RTMB::REPORT(uTot_g_nll)
  #   RTMB::ADREPORT(lambda_k)
  #
  #   print("test")
  #   ## ---------------------------------------------------------------
  #   ## 3. week × group counts  u_k_t  --------------------------------
  #   ## Build k_idx via Nsuper_for_joint
  #   ## Build observation matrix by joint group g × period t
  #   u_g_t <- array(0, dim = c(G,5,s))
  #   for (i in seq_along(n)) {
  #     t <- t_k[i]
  #     g <- sex[i]
  #     j <- t_l[i]
  #     u_g_t[g, j, t] <- u_g_t[g, j, t] + n[i]
  #   }
  #   RTMB::REPORT(u_g_t)
  #
  #   ## Multinomial likelihood over each joint group g
  #   uTot_g_t_nll <- 0
  #   for (g in 1:G) {
  #     for (j in 1:5) {
  #       prob <- (detect_mat[g,j, ] + 1e-6)
  #       prob <- prob / sum(prob)
  #       uTot_g_t_nll <- uTot_g_t_nll - dmultinom(
  #         x = u_g_t[g,j, ],
  #         prob = prob,
  #         log = TRUE
  #       )
  #     }
  #   }
  #   RTMB::REPORT(uTot_g_t_nll)
  #   RTMB::ADREPORT(prob)
  #
  #   nll <- sum(CJS_nll)
  #   nll <- nll + nll_u_w + nll_u_phi + nll_u_p + nll_u_t_vars
  #   nll <- nll + uTot_g_nll + uTot_g_t_nll
  #   #
  #   Ntotal <- sum(Nsuper_vec)
  #   RTMB::ADREPORT(Nsuper_vec)
  #   RTMB::ADREPORT(Ntotal)
  #
  #   return(nll)
  # }
  }
#
out <- list(rep = my_rep,
sd = my_sd)

saveRDS(out,"output/sex_loc_ttt_transition.rds")
#
for(i in 1:12){
    print(paste(i,round(sum(out$rep[[i]]$Nsuper_vec)), round(sum(out$rep[[i]]$Nsuper_vec - 1.96 * out$sd[[i]]$sd$Nsuper_vec)), round(sum(out$rep[[i]]$Nsuper_vec + 1.96 * out$sd[[i]]$sd$Nsuper_vec))))
  }
#
#
