model <- function(parms) {

  ## ------------------------------------------------------------------
  ## 0.  unpack data & parameters  ------------------------------------
  ## ------------------------------------------------------------------
  RTMB::getAll(data, parms)

  nll <- 0

  n_loc    <- if (!is.null(state)) length(unique(data[[state]])) else 1L
  RTMB::REPORT(n_loc)

  n_states <- 2 * n_loc + 1

  ## ------------------------------------------------------------------
  ## helpers -----------------------------------------------------------
  get_XZ_group <- function(proc, data, use_noperiod = FALSE) {
    print("proc")
    print(proc)
    # full‐size design matrices (n × p) and (n × q)
    X_full <- data$X_list[[proc]]
    Z_full <- data$Z_list[[proc]]
    print(X_full)
    print(Z_full)
    # choose which grouping factor to index
    gf <- if (use_noperiod) {
      data$design$group_noperiod[[proc]]
    } else {
      data$design$group_factor[[proc]]
    }
    # find the first row of each level
    idx <- match(levels(gf), gf)

    list(
      X = X_full[idx, , drop = FALSE],
      Z = Z_full[idx, , drop = FALSE]
    )
  }

  linpred <- function(g, xz, beta, u) {
    eta <- 0
    if (NCOL(xz$X)) eta <- eta + drop(xz$X[g, ] %*% beta)
    if (NCOL(xz$Z)) eta <- eta + drop(xz$Z[g, ] %*% u)
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

  K    <- nrow(xz$Nsuper$X)
  RTMB::REPORT(K)
  G <- nlevels(joint_group)
  print("G")
  print(G)
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
    nll_t_var <-  - sum(dnorm(u_t_var,0,1,TRUE))
  }
  RTMB::REPORT(M)
  RTMB::REPORT(u_t_var)

  ## helper ---------------------------------------------------------------
  taylor3_expm <- function(A) {
    I  <- diag(nrow(A))
    A2 <- A %*% A
    A3 <- A2 %*% A
    I + A + 0.5 * A2 + (1/6) * A3
  }

  ## ------------------------------------------------------------------
  ## 1. survival array  S (G × … × …) ---------------------------------
  S <- array(0, c(G, n_states, n_states))
  for (g in seq_len(G)) {
    # gphi <- phi_for_joint[g]                      # 1 … Gphi
    #
    # eta  <- linpred(gphi, xz$phi, beta_phi, u_phi)
    #
    # for (i in seq_len(n_loc)) {
    #   S[g, i, i]             <- -exp(eta)
    #   S[g, i, 2*n_loc + 1]   <-  exp(eta)
    # }
    # S[g,1:n_loc,1:n_loc] <- S[g,1:n_loc,1:n_loc] + M[1,,]
    # S[g,,] <- as.matrix(Matrix::expm(S[g,,]))
  }
  RTMB::REPORT(S)

  # --- Identify period in phi formula:
  print(design$formulas$phi)
  phi_has_period <- grepl(paste0("\\b","f_tk","\\b"),
                          paste(deparse(design$formulas$phi), collapse = ""))

  # --- Dimensions:
  s      <- max(t_k, na.omit(r_k))
  n_loc  <- if (!is.null(state)) length(levels(data[[state]])) else n_states
  n_states <- n_loc + 1L

  # --- Select the correct tiny-design for phi
  xz_phi <- get_XZ_group("phi", data, use_noperiod = !phi_has_period)
  G_use   <- nrow(xz_phi$X)  # equals nlevels(design$group_noperiod$phi) if no period,
  # else nlevels(design$group_factor$phi)

  # --- Precompute group→state lookup (for no-period or period case identical):
  grp_fac <- if (!phi_has_period) {
    design$group_noperiod[["phi"]]
  } else {
    design$group_factor[["phi"]]
  }
  # one representative row for each group level
  rep_idx <- match(levels(grp_fac), grp_fac)
  # pull the state values for each group
  if (!is.null(state)) {
    g_to_locs <- split(as.integer(data[[state]][rep_idx]),
                       seq_len(G_use))
  } else {
    # everyone maps to state 1
    g_to_locs <- lapply(seq_len(G_use), function(i) 1L)
  }

  # --- Allocate S2
  S2 <- array(0, dim = c(G_use, n_states, n_states, s))

  # --- Compute all φ‐rates in one shot
  linpred_vec <- function(g_idx, xz, beta, u) {
    eta <- numeric(length(g_idx))
    print(xz)
    if (ncol(xz$X) > 0) {
      eta <- eta + as.vector(xz$X[g_idx, , drop = FALSE] %*% beta)
    }
    if (ncol(xz$Z) > 0) {
      eta <- eta + as.vector(xz$Z[g_idx, , drop = FALSE] %*% u)
    }
    eta
  }
  eta_vec <- linpred_vec(seq_len(G_use), xz_phi, beta_phi, u_phi)

  if (!phi_has_period) {
    # ---- CASE 1: no f_tk in φ → replicate each Sg over all s
    for (g in seq_len(G_use)) {
      rate <- exp(eta_vec[g])

      # build one G_use×n_states generator
      Sg <- matrix(0, n_states, n_states)
      for (st in g_to_locs[[g]]) {
        Sg[st,     st    ] <- -rate
        Sg[st, n_states] <-  rate
      }
      Sg[1:n_loc, 1:n_loc] <- Sg[1:n_loc, 1:n_loc] + M[1, , ]
      Sg_exp <- as.matrix(Matrix::expm(Sg))

      # replicate across time
      for (tt in seq_len(s)) {
        S2[g, , , tt] <- Sg_exp
      }
    }

  } else {
    # ---- CASE 2: φ includes f_tk → fill only observed (g,tt) slices
    phi_full <- as.integer(design$group_factor[["phi"]])
    period_fac <- as.integer(data[[period]])

    # unique (g,tt) combos
    combos <- unique(
      data.frame(g = phi_full, tt = period_fac, stringsAsFactors = FALSE)
    )
    # for each combo build exactly one slice
    for (i in seq_len(nrow(combos))) {
      g  <- combos$g[i]
      tt <- combos$tt[i]
      rate <- exp(eta_vec[g])

      Sg <- matrix(0, n_states, n_states)
      for (st in g_to_locs[[g]]) {
        Sg[st,     st    ] <- -rate
        Sg[st, n_states] <-  rate
      }
      Sg[1:n_loc, 1:n_loc] <- Sg[1:n_loc, 1:n_loc] + M[1, , ]

      S2[g, , , tt] <- as.matrix(Matrix::expm(Sg))
    }
  }

  nll_u_phi <- 0
  if (NCOL(xz$phi$Z))                       # random-effect prior
    nll_u_phi <- - sum(dnorm(u_phi, 0, exp(log_sd_phi), TRUE))
  RTMB::REPORT(nll_u_phi)

  ## ------------------------------------------------------------------
  ## 2. detection array  D (G × … × …) --------------------------------
  design <- data$design
  fml    <- design$formulas$p
  has_tk <- grepl(paste0("\\b", 'f_tk', "\\b"), paste(deparse(fml), collapse=" "))

  # s        <- max(data$t_k, na.omit(data$r_k))
  n_loc    <- if (!is.null(state)) length(levels(data[[state]])) else data$n_states
  n_states <- n_loc + 1L

  # tiny‐design for p
  xz_p <- get_XZ_group("p", data, use_noperiod = !has_tk)
  G    <- nrow(xz_p$X)

  print(xz_p)

  # group→state lookup
  grp_fac  <- if (!has_tk) design$group_noperiod$p else design$group_factor$p
  rep_idx  <- match(levels(grp_fac), grp_fac)
  if (!is.null(state)) {
    loc_idx   <- as.integer(data[[state]][rep_idx])
    g_to_locs <- split(loc_idx, seq_len(G))
  } else {
    g_to_locs <- replicate(G, 1L, simplify = FALSE)
  }

  # compute all ηₚ’s and transform by logistic
  eta_p <- linpred_vec(seq_len(G), xz_p, beta_p, u_p)
  p_vec <- plogis(eta_p)

  # allocate
  print(c(G, n_states, n_states, data$s))
  D2 <- array(0, dim = c(G, n_states, n_states, s))
  if (!has_tk) {
    ## no‐period → same Dg every t
    for (g in seq_len(G)) {
      pg <- p_vec[g]
      Dg <- matrix(0, n_states, n_states)
      print(g_to_locs)
      print(g)
      print(seq_len(g_to_locs[[g]]))
      for (st in seq_len(g_to_locs[[g]])) {
        Dg[st,     st    ] <- 1 - pg   # remain undetected
        Dg[st, n_loc] <-     pg   # detected (moves to “observed” index)
      }
      Dg[n_states, n_states] <- 1      # dead stays dead
      for (tt in seq_len(s)) {
        D2[g,,,tt] <- Dg
      }
    }
  } else {
    ## period‐aware → only observed (g,tt)
    p_fac   <- as.integer(data$period_fac)
    combos  <- unique(data.frame(g = as.integer(grp_fac), tt = p_fac))
    for (i in seq_len(nrow(combos))) {
      g  <- combos$g[i]
      tt <- combos$tt[i]
      pg <- p_vec[g]

      D2[g,,,tt] <- diag(n_states)
      for (st in g_to_locs[[g]]) {
        D2[g, st,     st  , tt  ] <- 1 - pg
        D2[g, st, st + n_loc, tt] <-     pg
      }
    }
  }

  print(dim(D2))
  RTMB::REPORT(D)

#
#   ## ====================================================================
#   ## 2 bis — period-specific detection array  D2
#   ##         dim =  (G × n_states × n_states × s)
#   ## ====================================================================
#
#   D2 <- array(0, dim = c(G, n_states, n_states, s))
#
#   ## does the p-formula mention f_tk ?
#   p_has_period <- grepl("\\bf_tk\\b",
#                         paste(deparse(formulas$p), collapse = ""))
#
#   if (!p_has_period) {
#
#     ## ------------------------------------------------------------
#     ## case 1 :  p has NO period term  →  replicate ordinary D
#     ## ------------------------------------------------------------
#     for (tt in 1:s) D2[,,,tt] <- D      # D already built above
#
#   } else {
#
#     ## ------------------------------------------------------------
#     ## case 2 :  p DOES contain f_tk  →  build (g , tt) specific D
#     ## ------------------------------------------------------------
#
#     ## helper objects already available in your workspace -----------
#     period_fac   <- as.integer(data$period_fac)            # length n
#     p_grp_noper  <- design$group_noperiod[["p"]]           # length n
#
#     p_combo_df <- data.frame(
#       g    = as.integer(p_grp_noper),   # 1 … Gp  (without f_tk)
#       tt   = period_fac                # 1 … s
#     ) |>
#       dplyr::distinct()
#
#     RTMB::REPORT(p_combo_df)
#
#     for (k in seq_len(nrow(p_combo_df))) {
#
#       g  <- p_combo_df$g [k]      # group without period
#       tt <- p_combo_df$tt[k]      # period 1 … s
#
#       ## -------- linear predictor for this row --------------------
#       eta <- linpred(k, xz$p, beta_p, u_p)
#       p_i <- plogis(eta)
#
#       ## -------- tiny 2-state generator ---------------------------
#       D_tmp <- diag(n_states)
#       for (i in seq_len(n_loc)) {
#         D_tmp[i         , n_loc + i] <-  p_i         # live → “detected live”
#         D_tmp[i         , i        ] <-  1 - p_i     # live → live (not detected)
#       }
#
#       ## -------- store slice --------------------------------------
#       D2[g, , , tt] <- D_tmp
#     }
#   }
#   RTMB::REPORT(D2)

  nll_u_p <- 0
  if (NCOL(xz$p$Z))
    nll_u_p <-  - sum(dnorm(u_p, 0, exp(log_sd_p), TRUE))
  RTMB::REPORT(nll_u_p)
  RTMB::REPORT(u_p)

  ## ------------------------------------------------------------------
  ## 3. super-population predictor ------------------------------------
  eta_N <- numeric(K)
  for (k in seq_len(K))
    eta_N[k] <- linpred(k, xz$Nsuper, beta_Nsuper, u_Nsuper)
  Nsuper_vec <- exp(eta_N)
  RTMB::REPORT(eta_N)


  nll_u_Nsuper <- 0
  if (NCOL(xz$Nsuper$Z))
    nll_u_Nsuper <-  - sum(dnorm(u_Nsuper, 0, exp(log_sd_Nsuper), TRUE))
  RTMB::REPORT(nll_u_Nsuper)
  ## ------------------------------------------------------------------
  ##  temporal marking probability  v_w  (was missing)
  ## ------------------------------------------------------------------
  nll_u_v <- 0
  if(ncol(Z_list$v) == 0){
    v_w <- matrix(1 / s, nrow = G, ncol = s)
  } else{
    ## latent normals (one per primary period)
    nll_u_v <-  - sum(dnorm(u_v, 0, 1, log = TRUE))
    RTMB::REPORT(nll_u_v)
    ## transform:  Normal → Uniform → Beta
    v_u  <- pnorm(u_v)                               # length s
    # v_w  <- qbeta(v_u, shape1 = exp(v_a), shape2 = exp(v_b))
    v_w  <- qbeta(v_u, shape1 = exp(1), shape2 = exp(1))
  }
  RTMB::REPORT(v_w)


  has_fixed  <- NCOL(xz$w$X) > 0 && !all(colnames(xz$w$X) %in% "(Intercept)")
  has_random <- NCOL(xz$w$Z) > 0

  if (!has_fixed && !has_random) {
    # truly uniform model: w ~ 1
    pent_mat <- matrix(1/s, G_list$Gw, s)
  } else {
    # Build linear predictor
    Xw <- xz$w$X
    Zw <- xz$w$Z

    eta_fixed <- if (has_fixed)  Xw %*% beta_w else 0
    eta_random <- if (has_random) Zw %*% u_w else 0
    eta_w <- as.numeric(eta_fixed + eta_random)

    u_uniform <- pnorm(eta_w)
    w_group <- design$group_noperiod[["w"]]  # G_w groups
    periods <- as.integer(data$period_fac)   # s levels

    combo_df <- data.frame(
      group = as.integer(w_group),
      period = periods
    ) %>%
      distinct()
    RTMB::REPORT(combo_df)

    # create index for each valid (group, period) combo
    combo_df$row <- seq_len(nrow(combo_df))

    # Build a full matrix with NA initially
    delta_mat <- matrix(0, nrow = nlevels(w_group), ncol = s)

    # Get u_uniform and eta_w with same length as combo_df$row
    u_uniform  <- pnorm(eta_w)
    delta_vals <- qgamma(u_uniform, shape = exp(log_sd_w), scale = 1)

    for (i in seq_len(nrow(combo_df))) {
      g <- combo_df$group[i]
      t <- combo_df$period[i]
      delta_mat[g, t] <- delta_vals[i]
    }
    for(r in 1:nrow(delta_mat)){
      delta_mat[r,] <- delta_mat[r,]/sum((delta_mat[r,]))
    }
    RTMB::REPORT(delta_mat)

    ## map rows to joint groups
    w_levels <- levels(w_group)
    w_for_joint <- vapply(seq_len(G), function(g) {
      r <- which(joint_group == levels(joint_group)[g])[1]
      match(w_group[r], w_levels)
    }, integer(1))
    pent_mat <- delta_mat[w_for_joint, , drop = FALSE]


  }
  RTMB::REPORT(pent_mat)

  print("Done with S2")


  ## add prior
  nll_u_w <- 0
  if (has_random){
    nll_u_w <-  - sum(dnorm(u_w, 0, 1, log = TRUE))
  }
  RTMB::REPORT(nll_u_w)

  ## ------------------------------------------------------------------
  ## 6.  CJS likelihood (per-fish loop)  ------------------------------
  ## ------------------------------------------------------------------
  CJS_nll <- 0
  print(dim(S2))
  for (i in seq_along(t_k)) {

    entry      <- t_k[i]
    recapture  <- r_k[i]
    is_tagged  <- tag[i]
    count      <- n[i]

    ## probability a detected fish was tagged (temporal-marking)
    p_marked <- v_w[entry]

    ## seed α-vector
    alpha      <- numeric(n_states)
    idx_live   <- if (n_loc > 1) t_l[i] else 1L
    alpha[idx_live] <- 1

    ## ------- look-ups via new vectors -------------------------------
    g_joint <- joint_group[i]               # 1 … G  (store this in data)
    g_phi   <- phi_for_joint[g_joint]
    g_p     <- p_for_joint[g_joint]

    ## transition kernel for this fish
    Tmat <- S2[g_joint,,,t_k[i]] %*% D2[g_joint,,,t_k[i]]

    for (t in seq(entry , s - 1))
      alpha <- alpha %*% Tmat
      Tmat <- S2[g_joint,,,t+1] %*% D2[g_joint,,,t+1]


      # Set default in case of no recapture
      prob_rec <- 0

      # If recaptured, set detection index and extract from alpha
      if(is_tagged){
        if(n_loc>1){
          if (!is.na(r_l[i])) { #recapture
            idx_det <- if (n_loc > 1) n_loc + r_l[i] else n_loc + 1
            if (idx_det > length(alpha) || idx_det < 1)
              stop("idx_det out of bounds: ", idx_det)
            prob_rec <- alpha[idx_det]
          } else { #tagged but never recovered
            prob_rec <- sum(alpha[ (n_loc + 1):(2 * n_loc) ])
          }
        }else{
          prob_rec <- alpha[2]
        }
      }

        if (entry == s) prob_rec <- 0
    # print(paste(i, is_tagged, recapture))
    ## ------------------- likelihood contributions -------------------
    if (is_tagged) {
      CJS_nll <- CJS_nll - dbinom(count, count,  p_marked, log = TRUE)
      CJS_nll <- CJS_nll - dbinom( ifelse(is.na(recapture), 0, count),
                           count, prob_rec, log = TRUE)
    } else {
      # print("not tagged")
      CJS_nll <- CJS_nll - dbinom(count, count, 1 - p_marked, log = TRUE)
    }

  }
  RTMB::REPORT(CJS_nll)

  print("CJS")
  print(CJS_nll)
  ## ------------------------------------------------------------------
  ## 4.  calc_psi (unchanged – S,D,pent now already G-sized) ---------
  ## ------------------------------------------------------------------
  # print(G)
  calc_psi <- function(Smat, Dmat, pent_mat, sp_pr, s, n_loc) {
    G   <- nlevels(joint_group); n_states <- 2 * n_loc + 1
    psiPtot_g <- numeric(G); detect_mat <- matrix(0, G, s)

    if (n_loc == 1) sp_pr <- 1 else sp_pr <- sp_pr

    for (g in 1:G) {
      if(nlevels(w_for_joint)>1){
        init_t <- min(combo_df$period[combo_df$group == w_for_joint[g]])
      }else{
        init_t <- 1
      }

      alpha <- numeric(n_states); alpha[1:n_loc] <- pent_mat[w_for_joint[g],init_t] * sp_pr
      for (t in init_t:s) {
        ## entrants during week t
        for (j in 1:n_loc)
          alpha[j] <- alpha[j] +
            pent_mat[w_for_joint[g],t] * sp_pr[j] * (Smat[g,j,j,t] - 1) / log(Smat[g,j,j,t])
        p_vec   <- Dmat[g , 1:n_loc , 1:n_loc + n_loc,t]
        detect_mat[g,t] <- sum(alpha[1:n_loc] * p_vec)
        alpha <- alpha %*% Smat[g,,,t] %*% Dmat[g,,,t]
      }
      psiPtot_g[g] <- sum(detect_mat[g, ])
    }
    list(psiPtot_g=psiPtot_g, detect_mat=detect_mat)
  }

  out_g <- calc_psi(S2, D2, pent_mat, sp_pr, s, n_loc)
  psiPtot_g  <- out_g$psiPtot_g
  print(psiPtot_g)
  RTMB::REPORT(psiPtot_g)

  detect_mat <- out_g$detect_mat
  RTMB::REPORT(detect_mat)

  lambda_k <- Nsuper_vec[Nsuper_for_joint] * psiPtot_g           # length K
  RTMB::REPORT(lambda_k)
  RTMB::REPORT(Nsuper_vec)

  ## ---------------------------------------------------------------
  ## 2. group-wise total detections uTot_g -------------------------
  uTot_g <- data.frame(n = n, joint_group = joint_group) %>%
    group_by(joint_group) %>%
    summarise(uTot = sum(n), .groups = "drop")
  RTMB::REPORT(uTot_g)
  ## make sure it has length K:
  uTot_g_nll <- -sum(dpois(uTot_g$uTot, lambda_k, log = TRUE))
  RTMB::REPORT(uTot_g_nll)

  ## ---------------------------------------------------------------
  ## 3. week × group counts  u_k_t  --------------------------------
  ## Build k_idx via Nsuper_for_joint
  ## Build observation matrix by joint group g × period t
  u_g_t <- matrix(0, nrow = G, ncol = s)
  for (i in seq_along(n)) {
    g <- as.integer(joint_group[i])
    t <- t_k[i]
    u_g_t[g, t] <- u_g_t[g, t] + n[i]
  }
  RTMB::REPORT(u_g_t)


  pars <- list(beta_phi = beta_phi,
               beta_p = beta_p,
               u_t_var = u_t_var)

  RTMB::REPORT(pars)
  ## Multinomial likelihood over each joint group g
  print(detect_mat)
  uTot_g_t_nll <- 0
  for (g in seq_len(G)) {

    prob <- (detect_mat[g, ] + 1e-6)
    prob <- prob / sum(prob)
    uTot_g_t_nll <- uTot_g_t_nll - dmultinom(
      x = u_g_t[g, ],
      prob = prob,
      log = TRUE
    )
  }
  RTMB::REPORT(uTot_g_t_nll)

  nll <- CJS_nll
  nll <- nll + nll_u_w + nll_u_phi + nll_u_p + nll_u_v  + nll_t_var
  nll <- nll + uTot_g_nll + uTot_g_t_nll

  Ntotal <- sum(Nsuper_vec)
  RTMB::ADREPORT(Nsuper_vec)
  RTMB::ADREPORT(Ntotal)

  return(nll)
}
