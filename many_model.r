# rm(list=ls())
# library(dplyr)
# library(RTMB)
# library(tidyr)
# library(rlang)
#
#
#
# dans <- read.csv("data/NF Lewis_Chinook_abundance by year and run.csv")
# years <- 2013:2024
# my_rep <- list()
# my_sd <- list()
#
# f <- read.csv("C:/R_packages/kale/data/NF_Lewis_combined-2013_2024-2025-04-18.csv")

#ttt
# for(yyy in 1:1){
yyy <- 12
    d <- f %>%
      mutate(t_k = Period_Cap,
           r_k = Period_Recap,
           t_yr = Return_Yr) %>%
    filter(Return_Yr == years[yyy]) %>%
    mutate(t_l = if_else(is.na(TagState),NA,TagState),
           r_l = if_else(is.na(RecapState),NA,TagState)) %>%
    mutate(tag = ifelse(is.na(Tag1),FALSE,TRUE)) %>%
    mutate(f_yr = factor(t_yr, levels = unique(sort(as.numeric(t_yr))))) %>%
    mutate(f_tk = factor(t_k, levels = unique(sort(as.numeric(t_k))))) %>%
    mutate(f_tl = factor(t_l, levels = unique(sort(as.numeric(t_l))))) %>%
    mutate(f_rl = factor(r_l, levels = unique(sort(as.numeric(t_l))))) %>%
    mutate(f_FL = factor(FL, levels = unique(sort(as.numeric(FL))))) %>%
    mutate(f_sex = factor(Sex_final, levels = unique(Sex_final))) %>%
    group_by(t_k,r_k,t_l, r_l, tag, t_yr, f_yr, f_tk, f_rl, f_tl, f_sex) %>%
    summarise(n = n()) %>%
    filter((t_l<=r_l) | is.na(r_l))

  #Number of periods
  s <- max(d$t_k)

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
    r_l               = d$r_l,
    tag               = as.integer(d$tag),
    n                 = d$n,
    s                 = s,
    uTot_g            = uTot_g,
    sp_pr             = sp_pr
  )


  ttt_model <- function(parms) {

    ## ------------------------------------------------------------------
    ## 0.  unpack data & parameters  ------------------------------------
    ## ------------------------------------------------------------------
    RTMB::getAll(data, parms)

    nll <- 0

    n_loc    <- 1
    n_states <- 3
    ## ------------------------------------------------------------------
    ## 1. survival array  S (G × … × …) ---------------------------------
    S <- array(0, c(1, 3, 3, s))
    for (t in seq_len(s)) {
      S[1,1,1,t] <- -exp(phi + u_phi[t])
      S[1,1,3,t] <- -S[1,1,1,t]
      S[1,,,t] <- as.matrix(Matrix::expm(S[1,,,t]))
    }
    S_ <- S[1,1,1,]
    RTMB::REPORT(S)
    RTMB::ADREPORT(S_)
    nll_u_phi <- 0
    nll_u_phi <- -sum(dnorm(u_phi,0,exp(log_sd_phi),TRUE))

    D <- array(0, dim = c(1, 3, 3, s))
    for (t in 1:s) {
      D[1,1,1,t] <- -exp(p + u_p[t])
      D[1,1,2,t] <- -D[1,1,1,t]
      D[1,,,t] <- as.matrix(Matrix::expm(D[1,,,t]))
    }
    D_ <- D[1,1,1,]
    RTMB::REPORT(D)
    RTMB::ADREPORT(D_)
    nll_u_p <- 0
    nll_u_p <- -sum(dnorm(u_p,0,exp(log_sd_p),TRUE))

    Nsuper_vec <- exp(Nsuper)
    RTMB::REPORT(Nsuper_vec)

    delta_mat <- matrix(0, nrow = 1, ncol = s)

    # Get u_uniform and eta_w with same length as combo_df$row
    u_uniform  <- pnorm(u_w)
    delta_vals <- qgamma(u_uniform, shape = exp(log_sd_w), scale = 1)

    for (i in 1:s) {
      delta_mat[1, i] <- delta_vals[i]
    }
    delta_mat[1,] <- delta_mat[1,]/sum((delta_mat[1,]))
    RTMB::REPORT(delta_mat)

    pent_mat <- delta_mat[1, , drop = FALSE]
    RTMB::REPORT(pent_mat)
    RTMB::ADREPORT(pent_mat)

    ## add prior
    nll_u_w <- 0
    nll_u_w <-  - sum(dnorm(u_w, 0, 1, log = TRUE))
    RTMB::REPORT(nll_u_w)

    ## ------------------------------------------------------------------
    ## 6.  CJS likelihood (per-fish loop)  ------------------------------
    ## ------------------------------------------------------------------
    CJS_nll <- 0
    for (i in seq_along(t_k)) {

      entry      <- t_k[i]
      recapture  <- r_k[i]
      is_tagged  <- tag[i]
      count      <- n[i]


      ## seed α-vector
      alpha      <- numeric(3)
      idx_live   <- if (n_loc > 1) t_l[i] else 1L
      alpha[idx_live] <- 1

      ## transition kernel for this fish
      Tmat <- S[1,,,t_k[i]] %*% D[1,,,t_k[i]]
      if(entry<s){
        for (t in seq(entry , s - 1)){
          alpha <- alpha %*% Tmat
          # print(pste(i,t))
          Tmat <- S[1,,,t + 1] %*% D[1,,,t + 1]
        }
      }else{
        alpha <- alpha %*% Tmat
      }


      # Set default in case of no recapture
      prob_rec <- 0

      # If recaptured, set detection index and extract from alpha
      if(is_tagged){
        prob_rec <- alpha[2]
      }

      if (entry == s) prob_rec <- 0

      # print(paste(i, is_tagged, recapture))
      ## ------------------- likelihood contributions -------------------
      if (is_tagged) {
        CJS_nll <- CJS_nll - dbinom( ifelse(is.na(recapture), 0, count),
                                     count, prob_rec, log = TRUE)
      }
    }
    RTMB::REPORT(CJS_nll)
    ## ------------------------------------------------------------------
    ## 4.  calc_psi (unchanged – S,D,pent now already G-sized) ---------
    ## ------------------------------------------------------------------
    # print(G)
    calc_psi <- function(Smat, Dmat, pent_mat, sp_pr, s, n_loc) {
      G   <- 1; n_states <- 2 * n_loc + 1
      psiPtot_g <- numeric(G); detect_mat <- matrix(0, G, s)

      sp_pr <- 1
      init_t <- 1

      #Forward algorithm
      alpha <- numeric(n_states); alpha[1:n_loc] <- pent_mat[1,init_t] * sp_pr
      for (t in init_t:s) {
        ## entrants during week t
        for (j in 1:n_loc)
          alpha[j] <- alpha[j] +
            pent_mat[1,t] * sp_pr * (Smat[1,j,j,t] - 1) / log(Smat[1,j,j,t])
        p_vec   <- Dmat[1 , 1:n_loc , 1:n_loc + n_loc,t]
        detect_mat[1,t] <- sum(alpha[1:n_loc] * p_vec)
        alpha <- alpha %*% Smat[1,,,t] %*% Dmat[1,,,t]
      }
      psiPtot_g[1] <- sum(detect_mat[1, ])

      list(psiPtot_g=psiPtot_g, detect_mat=detect_mat)
    }

    out_g <- calc_psi(S, D, pent_mat, sp_pr, s, 1)
    psiPtot_g  <- out_g$psiPtot_g
    RTMB::REPORT(psiPtot_g)

    detect_mat <- out_g$detect_mat
    RTMB::REPORT(detect_mat)

    lambda_k <- Nsuper_vec * psiPtot_g           # length K
    RTMB::REPORT(lambda_k)
    RTMB::REPORT(Nsuper_vec)

    ## ---------------------------------------------------------------
    ## 2. group-wise total detections uTot_g -------------------------
    uTot_g <- data.frame(n = n, group = 1) %>%
      group_by(group) %>%
      summarise(uTot = sum(n), .groups = "drop")
    RTMB::REPORT(uTot_g)
    ## make sure it has length K:
    uTot_g_nll <- -sum(dpois(uTot_g$uTot, lambda_k, log = TRUE))
    RTMB::REPORT(uTot_g_nll)
    RTMB::ADREPORT(lambda_k)

    ## ---------------------------------------------------------------
    ## 3. week × group counts  u_k_t  --------------------------------
    ## Build k_idx via Nsuper_for_joint
    ## Build observation matrix by joint group g × period t
    u_g_t <- matrix(0, nrow = 1, ncol = s)
    for (i in seq_along(n)) {
      g <- 1
      t <- t_k[i]
      u_g_t[g, t] <- u_g_t[g, t] + n[i]
    }
    RTMB::REPORT(u_g_t)
    ## Multinomial likelihood over each joint group g
    uTot_g_t_nll <- 0
    for (g in 1:1) {

      prob <- (detect_mat[g, ] + 1e-6)
      prob <- prob / sum(prob)
      uTot_g_t_nll <- uTot_g_t_nll - dmultinom(
        x = u_g_t[g, ],
        prob = prob,
        log = TRUE
      )
    }
    RTMB::REPORT(uTot_g_t_nll)
    RTMB::ADREPORT(prob)

    nll <- CJS_nll
    nll <- nll + nll_u_w + nll_u_phi + nll_u_p
    nll <- nll + uTot_g_nll + uTot_g_t_nll

    Ntotal <- sum(Nsuper_vec)
    RTMB::ADREPORT(Nsuper_vec)
    RTMB::ADREPORT(Ntotal)

    return(nll)
  }

  params <- list(u_phi = rep(0,data$s),
                 log_sd_phi = -1,
                 u_p = rep(0,data$s),
                 log_sd_p = -1,
                 u_w = rep(0,data$s),
                 log_sd_w = -1,
                 phi = -1,
                 p = -1,
                 Nsuper = 9)

  ttt_obj <- MakeADFun(ttt_model,
                       data = data,
                       random = c("u_phi","u_p","u_w"),
                       parameters = params)

  opt <- nlminb(ttt_obj$par, ttt_obj$fn, ttt_obj$gr)

  my_rep[[yyy]] <- rep_ttt <- ttt_obj$report()
  sdr_ttt <- sdreport(ttt_obj)
  my_sd[[yyy]] <- list(est = as.list(sdr_ttt, "Estimate", report=TRUE),
                  sd = as.list(sdr_ttt, "Std. Error", report=TRUE))
# }
#
  out <- list(rep = my_rep,
              sd = my_sd)
# #
  saveRDS(out,"output/ttt.rds")
#
  for(i in 1:12){
    print(paste(i,round(out$rep[[i]]$Nsuper_vec), round(out$rep[[i]]$Nsuper_vec - 1.96 * out$sd[[i]]$sd$Nsuper_vec), round(out$rep[[i]]$Nsuper_vec + 1.96 * out$sd[[i]]$sd$Nsuper_vec), ",", dans$est[dans$Year%in%(2013:2024)][i]))
  }
#
#
