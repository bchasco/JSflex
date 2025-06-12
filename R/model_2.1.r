# model <- function(parms) {
#
#   ## ------------------------------------------------------------------
#   ## 0.  unpack data & parameters  ------------------------------------
#   ## ------------------------------------------------------------------
#   RTMB::getAll(data, parms)
#
#   nll <- 0
#
#   ## integer codes (1-based) – already length n
#   k_idx <- as.integer(Nsuper_group)   # 1 … K
#   g_idx <- as.integer(joint_group)    # 1 … G
#
#   # K <- nlevels(Nsuper_group)
#   # G <- nlevels(joint_group)
#
#   A <- Matrix::sparseMatrix(
#     i    = k_idx,
#     j    = g_idx,
#     x    = 1L,
#     dims = c(K, G)
#   )
#   A <- as.matrix(A)
#   A[A > 0] <- 1                       # turn counts back to 1’s
#   A <- as.matrix(A)
#   A_mat <- A
#
#   # phi_for_joint  <- phi_for_joint
#   # p_for_joint    <- p_for_joint
#   # w_for_joint    <- w_for_joint
#   # G              <- G                       # rows for S, D, pent
#
#   ## original objects (unchanged) ----------
#   s        <- s
#   n_loc    <- if (!is.null(state)) length(unique(data[[state]])) else 1L
#   n_states <- 2 * n_loc + 1
#
#   icnt <- 1
#   M <- array(0,c(1,n_loc,n_loc))
#   if(!is.null(state)){
#     for (i in seq_len(n_loc-1)) {
#       for (j in i:n_loc) {
#         if (j > i) {
#           M[1,i, j] <- exp(-u_t_var[icnt])
#           icnt <- icnt + 1L
#         }
#       }
#       # diagonal
#       M[1,i, i] <- -sum(M[1,i, -i])
#     }
#     nll <- nll - sum(dnorm(u_t_var,0,1,TRUE))
#   }
#   print(M)
#
#   ## ------------------------------------------------------------------
#   ## helpers -----------------------------------------------------------
#   get_XZ_group <- function(proc, data) {
#     X  <- data$X_list [[proc]]           # n × p (may be 0-col)
#     Z  <- data$Z_list [[proc]]           # n × q
#     gf <- data$group_factor[[proc]]      # length-n factor
#
#     idx <- match(levels(gf), gf)         # first row of every level
#     list(X = X[idx, , drop = FALSE],
#          Z = Z[idx, , drop = FALSE])
#   }
#
#   linpred <- function(g, xz, beta, u) {
#     eta <- 0
#     if (NCOL(xz$X)) eta <- eta + drop(xz$X[g, ] %*% beta)
#     if (NCOL(xz$Z)) eta <- eta + drop(xz$Z[g, ] %*% u)
#     eta                       # scalar
#   }
#
#   ## ------------------------------------------------------------------
#   ## design blocks -----------------------------------------------------
#   xz <- list(
#     phi     = get_XZ_group("phi",     data),
#     p       = get_XZ_group("p",       data),
#     w       = get_XZ_group("w",       data),
#     Nsuper  = get_XZ_group("Nsuper",  data)
#   )
#
#   Gphi <- nrow(xz$phi$X)              # #φ–groups
#   Gp   <- nrow(xz$p$X)
#   Gw   <- nrow(xz$w$X)
#   Gtvar <- nrow(xz$w$X)
#   K    <- nrow(xz$Nsuper$X)
#
#   ## ------------------------------------------------------------------
#   ## 1. survival array  S (G × … × …) ---------------------------------
#   S <- array(0, c(G, n_states, n_states))
#   for (g in seq_len(G)) {
#     gphi <- phi_for_joint[g]                      # 1 … Gphi
#
#     eta  <- linpred(gphi, xz$phi, beta_phi, u_phi)
#
#     for (i in seq_len(n_loc)) {
#       S[g, i, i]             <- -exp(eta)
#       S[g, i, 2*n_loc + 1]   <-  exp(eta)
#     }
#     S[g,1:n_loc,1:n_loc] <- S[g,1:n_loc,1:n_loc] + M[1,,]
#     S[g,,] <- as.matrix(Matrix::expm(S[g,,]))
#   }
#   print(S[1,,])
#   if (NCOL(xz$phi$Z))                       # random-effect prior
#     nll <- nll - sum(dnorm(u_phi, 0, exp(log_sd_phi), TRUE))
#
#   ## ------------------------------------------------------------------
#   ## 2. detection array  D (G × … × …) --------------------------------
#   D <- array(0, c(G, n_states, n_states))
#   for (g in seq_len(G)) {
#     gp   <- p_for_joint[g]                  # 1 … Gp
#     eta  <- linpred(gp, xz$p, beta_p, u_p)
#     p_i  <- plogis(eta)
#
#     D[g,,] <- diag(n_states)
#     for (i in seq_len(n_loc)) {
#       D[g, i, n_loc + i] <-  p_i
#       D[g, i, i        ] <- 1 - p_i
#     }
#   }
#
#   if (NCOL(xz$p$Z))
#     nll <- nll - sum(dnorm(u_p, 0, exp(log_sd_p), TRUE))
#
#   ## ------------------------------------------------------------------
#   ## 3. super-population predictor ------------------------------------
#   eta_N <- numeric(K)
#   for (k in seq_len(K))
#     eta_N[k] <- linpred(k, xz$Nsuper, beta_Nsuper, u_Nsuper)
#   Nsuper_vec <- exp(eta_N)
#
#   if (NCOL(xz$Nsuper$Z))
#     nll <- nll - sum(dnorm(u_Nsuper, 0, exp(log_sd_Nsuper), TRUE))
#
#   ## ------------------------------------------------------------------
#   ##  temporal marking probability  v_w  (was missing)
#   ## ------------------------------------------------------------------
#   if(ncol(Z_list$v) == 0){
#     v_w <- matrix(1 / s, nrow = G, ncol = s)
#   } else{
#     ## latent normals (one per primary period)
#     nll <- nll - sum(dnorm(u_v, 0, 1, log = TRUE))
#
#     ## transform:  Normal → Uniform → Beta
#     v_u  <- pnorm(u_v)                               # length s
#     # print(v_u)
#     # v_w  <- qbeta(v_u, shape1 = exp(v_a), shape2 = exp(v_b))
#     v_w  <- qbeta(v_u, shape1 = exp(1), shape2 = exp(1))
#   }
#
#   ## ------------------------------------------------------------------
#   ## 3.  entry probability  pent_mat  ---------------------------------
#   ##      • works for any RE term that *contains* the period factor
#   ##      • allows “missing” (group , period) combinations
#   ## ------------------------------------------------------------------
#   if (ncol(Z_list$w) == 0) {           # formula was   w ~ 1
#     pent_mat <- matrix(1 / s, G, s)    # uniform entry
#   } else {                             # random effects present
#     delta_mat            <- matrix(0, Gw, s)
#     G  <- nlevels(joint_group)
#     w_gf    <- group_factor$w
#     w_levels <- levels(w_gf)
#
#     ## for each joint‐level g, find its w‐level:
#     w_for_joint <- vapply(seq_len(G), function(g) {
#       # 1) which rows have this joint‐level?
#       this_joint <- levels(joint_group)[g]
#       r <- which(joint_group == this_joint)[1]
#       # 2) extract that row’s w‐group, then turn into an integer index
#       match(w_gf[r], w_levels)
#     }, integer(1))
#
#     if (anyNA(w_for_joint)) {
#       stop("Failed to map some joint‐levels back to a w‐level")
#     }
#
#
#     ## now pull the right rows out of delta_mat
#     pent_mat <- delta_mat[w_for_joint, , drop = FALSE]
#
#     ## --------------------------------------------------------------
#     ## 3·5  Gaussian-copula prior for the latent normals
#     ## --------------------------------------------------------------
#     nll <- nll - sum(dnorm(u_w, 0, 1, log = TRUE))
#   }
#   ## ------------------------------------------------------------------
#   ## 6.  CJS likelihood (per-fish loop)  ------------------------------
#   ## ------------------------------------------------------------------
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
#     Tmat <- S[g_joint,,] %*% D[g_joint,,]
#
#     for (t in seq(entry + 1, s))
#       alpha <- alpha %*% Tmat
#
#     ## probability of re-observation ---------------------------------
#     # print(alpha)
#     prob_rec <- 0
#     if (is_tagged) {
#       if (!is.na(recapture)) {              # recaptured
#         idx_det <- if (n_loc > 1) n_loc + r_l[i] else 2
#         prob_rec <- alpha[idx_det]
#       } else {                              # tagged but never recaptured
#         prob_rec <- sum(alpha[ (n_loc + 1):(2 * n_loc) ])
#       }
#     }
#     if (entry == s) prob_rec <- 0
#
#     # print(prob_rec)
#     ## ------------------- likelihood contributions -------------------
#     if (is_tagged) {
#       # print("count")
#       # print(count)
#       # print(p_marked)
#       nll <- nll - dbinom(count, count,  p_marked,     log = TRUE)
#       nll <- nll - dbinom( ifelse(is.na(recapture), 0, count),
#                            count, prob_rec, log = TRUE)
#     } else {
#       nll <- nll - dbinom(count, count, 1 - p_marked, log = TRUE)
#     }
#     # print(nll)
#   }
#
#   # print(nll)
#   ## ------------------------------------------------------------------
#   ## 4.  calc_psi (unchanged – S,D,pent now already G-sized) ---------
#   ## ------------------------------------------------------------------
#   # print(G)
#   calc_psi <- function(S, D, pent_mat, sp_pr, s, n_loc) {
#     G   <- dim(S)[1]; n_states <- 2 * n_loc + 1
#     psiPtot_g <- numeric(G); detect_mat <- matrix(0, G, s)
#
#     print(psiPtot_g)
#     if (n_loc == 1) sp_pr <- 1 else sp_pr <- sp_pr
#
#     for (g in 1:G) {
#       alpha <- numeric(n_states); alpha[1:n_loc] <- pent_mat[g,1] * sp_pr
#       for (t in 1:s) {
#         ## entrants during week t
#         for (j in 1:n_loc)
#           alpha[j] <- alpha[j] +
#             pent_mat[g,t] * sp_pr[j] * (S[g,j,j] - 1) / log(S[g,j,j])
#         p_vec   <- D[g , 1:n_loc , 1:n_loc + n_loc]
#         detect_mat[g,t] <- sum(alpha[1:n_loc] * p_vec)
#         alpha <- alpha %*% S[g,,] %*% D[g,,]
#       }
#       psiPtot_g[g] <- sum(detect_mat[g, ])
#     }
#     list(psiPtot_g=psiPtot_g, detect_mat=detect_mat)
#   }
#
#   out_g <- calc_psi(S, D, pent_mat, sp_pr, s, n_loc)
#   psiPtot_g  <- out_g$psiPtot_g
#   detect_mat <- out_g$detect_mat
#   ## ------------------------------------------------------------------
#   ## 5.  aggregate to Nsuper via pre-built A -------------------------
#   ## ------------------------------------------------------------------
#   psi_dag <- A_mat %*% psiPtot_g                 # K × 1
#
#
#   lambda_k <- Nsuper_vec * psi_dag           # length K
#
#   ## ---------------------------------------------------------------
#   ## 2. group-wise total detections uTot_g -------------------------
#   uTot_g <- as.numeric(tapply(n, k_idx, sum))
#   ## make sure it has length K:
#   if (length(uTot_g) < K) uTot_g <- c(uTot_g, rep(0, K - length(uTot_g)))
#   nll <- nll - sum(dpois(uTot_g, lambda_k, log = TRUE))
#
#   ## ---------------------------------------------------------------
#   ## 3. week × group counts  u_k_t  --------------------------------
#   u_k_t <- matrix(0, K, s)
#   u_k_t[cbind(k_idx, t_k)] <- n + u_k_t[cbind(k_idx, t_k)]
#
#   tmp  <- matrix(0, K, s)
#   for (i in seq_along(t_k)) {
#     kk <- k_idx[i]
#     tt <- t_k  [i]
#     tmp[kk, tt] <- tmp[kk, tt] + n[i]
#   }
#   u_k_t2 <- tmp
#
#   ## detection by week & group (already built)
#   detect_k_t <- A %*% detect_mat            # K × s
#   multP_k_t  <- detect_k_t / rowSums(detect_k_t)
#
#   for (k in 1:K)
#     nll <- nll - dmultinom(u_k_t2[k, ]/sum(u_k_t2[k, ]), prob = multP_k_t[k, ], log = TRUE)
#
#   RTMB::REPORT(S)
#   RTMB::REPORT(D)
#   RTMB::REPORT(Gw)
#   RTMB::REPORT(Gphi)
#   RTMB::REPORT(Gp)
#   RTMB::REPORT(G)
#   RTMB::REPORT(K)
#   RTMB::REPORT(s)
#   RTMB::REPORT(pent_mat)
#   RTMB::REPORT(psiPtot_g)
#   RTMB::REPORT(psi_dag)
#   RTMB::REPORT(detect_mat)
#   RTMB::REPORT(Nsuper_vec)
#   RTMB::REPORT(uTot_g)
#   RTMB::REPORT(lambda_k)
#   RTMB::REPORT(detect_k_t)
#   RTMB::REPORT(A_mat)
#   RTMB::REPORT(A)
#   RTMB::REPORT(xz)
#   RTMB::REPORT(u_k_t)
#   RTMB::REPORT(multP_k_t)
#   RTMB::REPORT(detect_k_t)
#   RTMB::REPORT(u_k_t2)
#   RTMB::REPORT(M)
#   RTMB::REPORT(n_loc)
#   # psiPtot_g=psiPtot_g, psi_dag=psi_dag))
#   return(nll)
# }
