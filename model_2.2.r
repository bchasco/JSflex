model <- function(parms) {

  ## ------------------------------------------------------------------
  ## 0.  unpack data & parameters  ------------------------------------
  ## ------------------------------------------------------------------
  RTMB::getAll(data,parms)

  nll <- 0

  n_loc    <- if (!is.null(state)) length(unique(data[[state]])) else 1L
  RTMB::REPORT(n_loc)

  n_states <- 2 * n_loc + 1
  # --- Dimensions:
  n_time      <- max(t_k, na.omit(r_k))

  #
  ## ------------------------------------------------------------------
  get_XZ_group <- function(proc, data) {
    if (!design$has_period[[proc]]) {
      # — pull from the *obs* matrices and then reduce them —
      X_obs <- design$X_obs[[proc]]
      Z_obs <- design$Z_obs[[proc]]
      gf    <- design$group_noperiod[[proc]]
      # find one representative row per level of the no-period grouping:
      print(nlevels(gf))
      print(gf)
      rep_idx <- match(levels(gf), gf)
      X_min <- X_obs[rep_idx, , drop = FALSE]
      Z_min <- if (ncol(Z_obs)>0) Z_obs[rep_idx, , drop = FALSE] else Z_obs

    } else {
      # — full grouping: X_group/Z_group are *already* minimal —
      X_min <- data$design$X_group[[proc]]
      Z_min <- data$design$Z_group[[proc]]
    }

    list(X = X_min, Z = Z_min)
  }

  linpred <- function(xz, beta, u) {
    eta <- 0
    if (NCOL(xz$X)) eta <- eta + drop(xz$X %*% beta)
    if (NCOL(xz$Z)) eta <- eta + drop(xz$Z %*% u)
    eta                       # scalar
  }

  ## helper
  linpred_period <- function(g, tt, X, Z,
                                 phi_group, period_fac,
                                 beta, u) {
    row_id <- which(phi_group == g & period_fac == tt)[1]
    eta <- 0
    if (NCOL(X)) eta <- eta + drop(X[row_id, ] %*% beta)
    if (NCOL(Z)) eta <- eta + drop(Z[row_id, ] %*% u)
    eta
  }

  ## ------------------------------------------------------------------
  ## design blocks -----------------------------------------------------
  xz <- list(
    phi     = get_XZ_group("phi",     data),
    p       = get_XZ_group("p",       data),
    w       = get_XZ_group("w",       data),
    t_var       = get_XZ_group("t_var",       data),
    Nsuper  = get_XZ_group("Nsuper",  data)
  )
  RTMB::REPORT(xz)

  #
  G_list <- list(
    Gphi   = nlevels(design$group_noperiod[["phi"]]),
    Gp     = nlevels(design$group_noperiod[["p"]]),
    Gw     = nlevels(design$group_noperiod[["w"]]),
    Gtvar  = nlevels(design$group_noperiod[["t_var"]])
  )
  RTMB::REPORT(G_list)
  #
  K    <- nrow(xz$Nsuper$X)
  RTMB::REPORT(K)

  G <- nlevels(joint_group)
  RTMB::REPORT(G)

  icnt <- 1
  nll_t_var <- 0
  M <- array(0,c(1,n_loc,n_loc))
  if(!is.null(state)){
    for (i in seq_len(n_loc-1)) {
      for (j in i:n_loc) {
        if (j > i) {
          M[1,i, j] <- exp(-u_t_var[icnt])
          icnt <- icnt + 1L
        }
      }
      # diagonal
      M[1,i, i] <- -sum(M[1,i, -i])
    }
    nll_t_var <-  - sum(dnorm(u_t_var,0, u_sd_t_var,TRUE))
  }
  RTMB::REPORT(M)
  RTMB::REPORT(u_t_var)


  # --- Select the correct tiny-design for phi
  xz_phi <- get_XZ_group("phi", data)
  eta_phi <- linpred(xz$phi, beta_phi, u_phi)
  # --- Allocate S2
  S2 <- array(0, dim = c(length(unique(lk$phi$lookup$g)), n_states, n_states, n_time))
  for(i in 1:nrow(lk$phi$lookup)){
    g <- lk$phi$lookup$g[i]       # 1st dim
    s <- lk$phi$lookup$state[i]   # 2nd dim
    t <- lk$phi$lookup$time[i]   # 3rd dim (survival keeps you in the same state)
    eta_idx <- lk$phi$lookup$design_row[i]
    S2[g,s,s,t] <- eta_phi[eta_idx]
  }
  RTMB::REPORT(S2)
  RTMB::REPORT(eta_phi)


  # --- Select the correct tiny-design for phi
  xz_p <- get_XZ_group("p", data)
  eta_p <- linpred(xz$p, beta_p, u_p)
  # --- Allocate S2
  D2 <- array(0, dim = c(length(unique(lk$p$lookup$g)), n_states, n_states, n_time))
  for(i in 1:nrow(lk$phi$lookup)){
    g <- lk$p$lookup$g[i]       # 1st dim
    s <- lk$p$lookup$state[i]   # 2nd dim
    t <- lk$p$lookup$time[i]   # 3rd dim (survival keeps you in the same state)
    eta_idx <- lk$p$lookup$design_row[i]
    D2[g,s,s,t] <- eta_p[eta_idx]
    if(!design$has_period[["p"]]){
      D2[g,s,s,] <- D2[g,s,s,1]
    }
  }
  RTMB::REPORT(D2)
  RTMB::REPORT(eta_p)


  # --- Select the correct tiny-design for phi
  xz_Nsuper <- get_XZ_group("Nsuper", data)
  Nsuper <- array(0, dim = c(length(unique(lk$Nsuper$lookup$g)), n_states, n_states, n_time))
  RTMB::REPORT(xz_Nsuper)
  RTMB::REPORT(Nsuper)

  xz_w <- get_XZ_group("w", data)
  w <- array(0, dim = c(length(unique(lk$w$lookup$g)), n_states, n_states, n_time))
  RTMB::REPORT(xz_w)
  RTMB::REPORT(w)

  xz_t_var <- get_XZ_group("t_var", data)
  t_var <- array(0, dim = c(length(unique(lk$t_var$lookup$g)), n_states, n_states, n_time))
  RTMB::REPORT(xz_t_var)
  RTMB::REPORT(t_var)

  ## ------------------------------------------------------------------
  ## 3. super-population predictor ------------------------------------
  # eta_N <- numeric(K)
  # for (k in seq_len(K))
  #   eta_N[k] <- linpred(k, xz$Nsuper, beta_Nsuper, u_Nsuper)
  # Nsuper_vec <- exp(eta_N)
  # RTMB::REPORT(eta_N)
#
#
#   nll_u_Nsuper <- 0
#   if (NCOL(xz$Nsuper$Z))
#     nll_u_Nsuper <-  - sum(dnorm(u_Nsuper, 0, exp(log_sd_Nsuper), TRUE))
#   RTMB::REPORT(nll_u_Nsuper)
#   ## ------------------------------------------------------------------
#   ##  temporal marking probability  v_w  (was missing)
#   ## ------------------------------------------------------------------
#   nll_u_v <- 0
#   if(ncol(Z_list$v) == 0){
#     v_w <- matrix(1 / s, nrow = G, ncol = s)
#   } else{
#     ## latent normals (one per primary period)
#     nll_u_v <-  - sum(dnorm(u_v, 0, 1, log = TRUE))
#     RTMB::REPORT(nll_u_v)
#     ## transform:  Normal → Uniform → Beta
#     v_u  <- pnorm(u_v)                               # length s
#     # v_w  <- qbeta(v_u, shape1 = exp(v_a), shape2 = exp(v_b))
#     v_w  <- qbeta(v_u, shape1 = exp(1), shape2 = exp(1))
#   }
#   RTMB::REPORT(v_w)
#
#
#   has_fixed  <- NCOL(xz$w$X) > 0 && !all(colnames(xz$w$X) %in% "(Intercept)")
#   has_random <- NCOL(xz$w$Z) > 0
#
#   if (!has_fixed && !has_random) {
#     # truly uniform model: w ~ 1
#     pent_mat <- matrix(1/s, G_list$Gw, s)
#   } else {
#     # Build linear predictor
#     Xw <- xz$w$X
#     Zw <- xz$w$Z
#
#     eta_fixed <- if (has_fixed)  Xw %*% beta_w else 0
#     eta_random <- if (has_random) Zw %*% u_w else 0
#     eta_w <- as.numeric(eta_fixed + eta_random)
#
#     u_uniform <- pnorm(eta_w)
#     w_group <- design$group_noperiod[["w"]]  # G_w groups
#     periods <- as.integer(data$period_fac)   # s levels
#
#     combo_df <- data.frame(
#       group = as.integer(w_group),
#       period = periods
#     ) %>%
#       distinct()
#     RTMB::REPORT(combo_df)
#
#     # create index for each valid (group, period) combo
#     combo_df$row <- seq_len(nrow(combo_df))
#
#     # Build a full matrix with NA initially
#     delta_mat <- matrix(0, nrow = nlevels(w_group), ncol = s)
#
#     # Get u_uniform and eta_w with same length as combo_df$row
#     u_uniform  <- pnorm(eta_w)
#     delta_vals <- qgamma(u_uniform, shape = exp(log_sd_w), scale = 1)
#
#     for (i in seq_len(nrow(combo_df))) {
#       g <- combo_df$group[i]
#       t <- combo_df$period[i]
#       delta_mat[g, t] <- delta_vals[i]
#     }
#     for(r in 1:nrow(delta_mat)){
#       delta_mat[r,] <- delta_mat[r,]/sum((delta_mat[r,]))
#     }
#     RTMB::REPORT(delta_mat)
#
#     ## map rows to joint groups
#     w_levels <- levels(w_group)
#     w_for_joint <- vapply(seq_len(G), function(g) {
#       r <- which(joint_group == levels(joint_group)[g])[1]
#       match(w_group[r], w_levels)
#     }, integer(1))
#     pent_mat <- delta_mat[w_for_joint, , drop = FALSE]
#
#
#   }
#   RTMB::REPORT(pent_mat)
#
#   print("Done with S2")
#
#
#   ## add prior
#   nll_u_w <- 0
#   if (has_random){
#     nll_u_w <-  - sum(dnorm(u_w, 0, 1, log = TRUE))
#   }
#   RTMB::REPORT(nll_u_w)
#
#   ## ------------------------------------------------------------------
#   ## 6.  CJS likelihood (per-fish loop)  ------------------------------
#   ## ------------------------------------------------------------------
#   CJS_nll <- 0
#   print(dim(S2))
#   for (i in seq_along(t_k)) {
#
#     entry      <- t_k[i]
#     recapture  <- r_k[i]
#     is_tagged  <- tag[i]
#     count      <- n[i]
#
#     ## probability a detected fish was tagged (temporal-marking)
#     p_marked <- v_w[entry]
#
#     ## seed α-vector
#     alpha      <- numeric(n_states)
#     idx_live   <- if (n_loc > 1) t_l[i] else 1L
#     alpha[idx_live] <- 1
#
#     ## ------- look-ups via new vectors -------------------------------
#     g_joint <- joint_group[i]               # 1 … G  (store this in data)
#     g_phi   <- phi_for_joint[g_joint]
#     g_p     <- p_for_joint[g_joint]
#
#     ## transition kernel for this fish
#     Tmat <- S2[g_joint,,,t_k[i]] %*% D2[g_joint,,,t_k[i]]
#
#     for (t in seq(entry , s - 1))
#       alpha <- alpha %*% Tmat
#       Tmat <- S2[g_joint,,,t+1] %*% D2[g_joint,,,t+1]
#
#
#       # Set default in case of no recapture
#       prob_rec <- 0
#
#       # If recaptured, set detection index and extract from alpha
#       if(is_tagged){
#         if(n_loc>1){
#           if (!is.na(r_l[i])) { #recapture
#             idx_det <- if (n_loc > 1) n_loc + r_l[i] else n_loc + 1
#             if (idx_det > length(alpha) || idx_det < 1)
#               stop("idx_det out of bounds: ", idx_det)
#             prob_rec <- alpha[idx_det]
#           } else { #tagged but never recovered
#             prob_rec <- sum(alpha[ (n_loc + 1):(2 * n_loc) ])
#           }
#         }else{
#           prob_rec <- alpha[2]
#         }
#       }
#
#         if (entry == s) prob_rec <- 0
#     # print(paste(i, is_tagged, recapture))
#     ## ------------------- likelihood contributions -------------------
#     if (is_tagged) {
#       CJS_nll <- CJS_nll - dbinom(count, count,  p_marked, log = TRUE)
#       CJS_nll <- CJS_nll - dbinom( ifelse(is.na(recapture), 0, count),
#                            count, prob_rec, log = TRUE)
#     } else {
#       # print("not tagged")
#       CJS_nll <- CJS_nll - dbinom(count, count, 1 - p_marked, log = TRUE)
#     }
#
#   }
#   RTMB::REPORT(CJS_nll)
#
#   print("CJS")
#   print(CJS_nll)
#   ## ------------------------------------------------------------------
#   ## 4.  calc_psi (unchanged – S,D,pent now already G-sized) ---------
#   ## ------------------------------------------------------------------
#   # print(G)
#   calc_psi <- function(Smat, Dmat, pent_mat, sp_pr, s, n_loc) {
#     G   <- nlevels(joint_group); n_states <- 2 * n_loc + 1
#     psiPtot_g <- numeric(G); detect_mat <- matrix(0, G, s)
#
#     if (n_loc == 1) sp_pr <- 1 else sp_pr <- sp_pr
#
#     for (g in 1:G) {
#       if(nlevels(w_for_joint)>1){
#         init_t <- min(combo_df$period[combo_df$group == w_for_joint[g]])
#       }else{
#         init_t <- 1
#       }
#
#       alpha <- numeric(n_states); alpha[1:n_loc] <- pent_mat[w_for_joint[g],init_t] * sp_pr
#       for (t in init_t:s) {
#         ## entrants during week t
#         for (j in 1:n_loc)
#           alpha[j] <- alpha[j] +
#             pent_mat[w_for_joint[g],t] * sp_pr[j] * (Smat[g,j,j,t] - 1) / log(Smat[g,j,j,t])
#         p_vec   <- Dmat[g , 1:n_loc , 1:n_loc + n_loc,t]
#         detect_mat[g,t] <- sum(alpha[1:n_loc] * p_vec)
#         alpha <- alpha %*% Smat[g,,,t] %*% Dmat[g,,,t]
#       }
#       psiPtot_g[g] <- sum(detect_mat[g, ])
#     }
#     list(psiPtot_g=psiPtot_g, detect_mat=detect_mat)
#   }
#
#   out_g <- calc_psi(S2, D2, pent_mat, sp_pr, s, n_loc)
#   psiPtot_g  <- out_g$psiPtot_g
#   print(psiPtot_g)
#   RTMB::REPORT(psiPtot_g)
#
#   detect_mat <- out_g$detect_mat
#   RTMB::REPORT(detect_mat)
#
#   lambda_k <- Nsuper_vec[Nsuper_for_joint] * psiPtot_g           # length K
#   RTMB::REPORT(lambda_k)
#   RTMB::REPORT(Nsuper_vec)
#
#   ## ---------------------------------------------------------------
#   ## 2. group-wise total detections uTot_g -------------------------
#   uTot_g <- data.frame(n = n, joint_group = joint_group) %>%
#     group_by(joint_group) %>%
#     summarise(uTot = sum(n), .groups = "drop")
#   RTMB::REPORT(uTot_g)
#   ## make sure it has length K:
#   uTot_g_nll <- -sum(dpois(uTot_g$uTot, lambda_k, log = TRUE))
#   RTMB::REPORT(uTot_g_nll)
#
#   ## ---------------------------------------------------------------
#   ## 3. week × group counts  u_k_t  --------------------------------
#   ## Build k_idx via Nsuper_for_joint
#   ## Build observation matrix by joint group g × period t
#   u_g_t <- matrix(0, nrow = G, ncol = s)
#   for (i in seq_along(n)) {
#     g <- as.integer(joint_group[i])
#     t <- t_k[i]
#     u_g_t[g, t] <- u_g_t[g, t] + n[i]
#   }
#   RTMB::REPORT(u_g_t)
#
#
#   pars <- list(beta_phi = beta_phi,
#                beta_p = beta_p,
#                u_t_var = u_t_var)
#
#   RTMB::REPORT(pars)
#   ## Multinomial likelihood over each joint group g
#   print(detect_mat)
#   uTot_g_t_nll <- 0
#   for (g in seq_len(G)) {
#
#     prob <- (detect_mat[g, ] + 1e-6)
#     prob <- prob / sum(prob)
#     uTot_g_t_nll <- uTot_g_t_nll - dmultinom(
#       x = u_g_t[g, ],
#       prob = prob,
#       log = TRUE
#     )
#   }
#   RTMB::REPORT(uTot_g_t_nll)
#
#   nll <- CJS_nll
#   nll <- nll + nll_u_w + nll_u_phi + nll_u_p + nll_u_v  + nll_t_var
#   nll <- nll + uTot_g_nll + uTot_g_t_nll
#
#   Ntotal <- sum(Nsuper_vec)
#   RTMB::ADREPORT(Nsuper_vec)
#   RTMB::ADREPORT(Ntotal)

  return(0)
}
