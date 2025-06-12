# library(RTMB)
# model <- function(parms){
#
#   RTMB::getAll(data,parms)
#   ## Transformed parameters
#
#   Ntot <- exp(beta_Nsuper)
#
#   ## Negative log-likelihood
#   nll <- 0
#
#   out <- list()
#
#   #Temporal marking probability
#   nll <- nll - sum(dnorm(logit_v, 0, 1, log = TRUE));       # Assign N(0,1) distribution u
#   v_u <-  pnorm(logit_v, 0, 1);  # Uniformly distributed variables (on [0,1])
#   out$v_w <- v_w <- qbeta(v_u, shape1 = exp(v_a), shape2 = exp(v_b));
#
#   # For other components
#   n_loc <- if (!is.null(state)) length(unique(data[[state]])) else 1L
#
#   map_groups <- function(from_group, to_group) {
#     stopifnot(length(from_group) == length(to_group))
#     lookup <- data.frame(from = from_group, to = to_group)
#     names(lookup) <- c("from", "to")
#     map <- aggregate(to ~ from, data = lookup, FUN = function(x) {
#       ux <- unique(x)
#       if (length(ux) == 1) return(ux)
#       ## fallback to mode if inconsistent mapping
#       ux[which.max(tabulate(match(x, ux)))]
#     })
#
#     out <- integer(max(map$from))
#     out[map$from] <- map$to
#     return(out)
#   }
#
#   match_group_by_value <- function(ref_group, target_groups) {
#     for (name in names(target_groups)) {
#       tgt <- target_groups[[name]]
#       if (length(ref_group) != length(tgt)) next
#       # Match up to relabeling
#       if (all(rank(ref_group, ties.method = "first") == rank(tgt, ties.method = "first"))) {
#         return(name)
#       }
#     }
#     return(NULL)
#   }
#
#   #All the intermediate variables
#   if(!is.null(struct$groups$v$group_v_re1)){
#     ng_v <- max(struct$groups$v$group_v_re1)
#   }else{
#     ng_v <- 1
#   }
#
#   lambda <- rep(0,s)
#   temp <- rep(0,s)
#   tau <- rep(0,s)
#   psi <- rep(0,s)
#   multP <- rep(0,s)
#
#   n_states <- n_loc + n_loc + 1  # live + detected + dead
#   M_group   <- struct$groups$t_var
#   phi_group <- struct$groups$phi
#   p_group   <- struct$groups$p
#   w_group   <- struct$groups$w
#
#   # Movement matrix: uniform probability of moving to any equal or higher location
#   icnt <- 1
#   if(!is.null(struct$groups$t_var$group_t_var_re1)){
#     ng_M <- max(struct$groups$t_var$group_t_var_re1)
#   }else{
#     ng_M <- 1
#   }
#   M <- array(0, c(ng_M, n_loc, n_loc))
#   print(t_var)
#   if(n_loc > 1){
#     nll <- nll - sum(dnorm(t_var,0,1,TRUE))
#     out$t_var <- t_var
#     for(g in 1:ng_M){
#       for (i in 1:(n_loc-1)) {
#         for (j in i:n_loc) {
#           if(j>i){
#             M[g, i, j] <- exp(-(t_var[icnt]))# * exp(M_sigma)
#             icnt <- icnt + 1
#           }
#         }
#         M[g,i, i] <- -sum(M[g,i,])
#       }
#     }
#   }
#   out$M <- M
#   out$M_for_phi_group <- M_for_phi_group  <- map_groups(phi_group, M_group)
#   print(M_for_phi_group)
#   print(dim(M))
#   if(!is.null(struct$groups$phi$group_phi_re1)){
#     ng_phi <- max(struct$groups$phi$group_phi_re1)
#   }else{
#     ng_phi <- 1
#   }
#
#   S <- array(0, c(ng_phi, n_states, n_states))
#   for(g in 1:ng_phi){
#     for(i in 1:n_loc){
#       if(!is.null(struct$groups$phi$group_phi_re1)){
#         eta <- beta_phi[1] + u_phi[g]
#       }else{
#         eta <- beta_phi[1]
#       }
#
#       # print(g)
#       # print(dim(S[g,,]))
#       S[g, i, i] <- -exp(eta)
#       S[g, i, 2 * n_loc + 1] <- exp(eta)
#     }
#
#     if(n_loc>1){
#       S[g, 1:n_loc, 1:n_loc] <- S[g,1:n_loc, 1:n_loc] + M[M_for_phi_group[g],,]
#     }
#     # S[cbind(1:n_loc, 2 * n_loc + 1)] <- exp(phi)  # Death to undetected
#     S[g,,] <- as.matrix(Matrix::expm(S[g,,]))
#   }
#   out$S <- S
#
#   # Detection matrix
#   if(!is.null(struct$groups$p$group_p_re1)){
#     ng_p <- max(struct$groups$p$group_p_re1)
#   }else{
#     ng_p <- 1
#   }
#   ## Detection matrix (one row per p‑group)
#   D <- array(0, c(ng_p, n_states, n_states))
#   for (g in 1:ng_p) {
#     D[g,,] <- diag(n_states)        # <‑‑ keep absorbing rows!
#
#     for (i in 1:n_loc) {
#       eta <- if (ng_p > 1) beta_p[1] + u_p[g] else beta_p[1]
#       p_i <- plogis(eta)
#
#       D[g, i, n_loc + i] <-  p_i         # live → detected
#       D[g, i, i        ] <- 1 - p_i      # live stays undetected
#     }
#   }
#   out$D <- D
#
#   ## -------------------------------------------------------------
#   ## entry probability  (week × year specific)  ------------------
#   ## -------------------------------------------------------------
#   if (length(struct$Z_list$w) == 0L) {
#     ## formula was ~1 (uniform entry)  ----------------------------
#     pent_mat <- matrix(1 / s, nrow = 1, ncol = s)     # every row uniform
#   } else {
#
#     n_row_re <- length(u_w) %/% s          # ALWAYS works
#     # print(n_ro)
#     # stopifnot(n_row_re * s == length(u_w)) # sanity
#
#     ## reshape: rows = w-groups (year, state, …), cols = weeks
#     u_w_mat <- matrix(u_w, nrow = n_row_re, ncol = s, byrow = TRUE)
#
#     nll <- nll - sum(dnorm(u_w_mat, 0, 1, TRUE))
#     eta_mat     <- beta_w + u_w_mat
#     delta_u_mat <- pnorm(eta_mat)
#
#     delta_w_mat <- matrix(
#       qgamma(delta_u_mat,
#              shape = exp(beta_w), scale = 1),
#       nrow = n_row_re,
#       ncol = s,
#       byrow = TRUE)
#
#     if (n_row_re > 1) {
#       pent_mat <- delta_w_mat / rowSums(delta_w_mat)   # each row sums to 1
#     } else {
#       pent_mat <- delta_w_mat / sum(delta_w_mat)
#     }
#
#   }
#   out$pent_mat <- pent_mat
#   ## ----------  harmonise rows with φ-groups  --------------------
#   #pent only interacts with S in the allpha likelihood, so only the S, not the D, dim is important
#   G_phi <- dim(S)[1]                             # rows in survival array
#   if (nrow(pent_mat) < G_phi) {                  # replicate if needed
#     pent_mat <- pent_mat[rep(1, G_phi), , drop = FALSE]
#   }
#   if (nrow(pent_mat) > G_phi) {                  # or trim (rare case)
#     pent_mat <- pent_mat[seq_len(G_phi), , drop = FALSE]
#   }
#
#   for (i in 1:length(t_k)) {#individual
#     entry <- t_k[i]
#     recapture <- r_k[i]
#     is_tagged <- tag[i]
#     count <- n[i]
#
#     p_marked <- v_w[entry]
#     alpha <- matrix(0, nrow = 1, ncol = n_states)
#     # alpha[ t_l[i] ] <- 1  # starts alive at its tagging location
#     ## index of the individual’s live state at first detection
#     idx_live  <- if (n_loc > 1) t_l[i] else 1L
#     idx_group <- struct$groups$phi[[1]][i]   # φ–group row
#
#     alpha              <- numeric(n_states)
#     alpha[idx_live]    <- 1                  # seed
#
#     Tmat <- S[idx_group , , ] %*% D[struct$groups$p[[1]][i] , , ]
#
#     for (t in (t_k[i] + 1):s)
#       alpha <- alpha %*% Tmat
#
#     # Set default in case of no recapture
#     prob_rec <- 0
#
#     # If recaptured, set detection index and extract from alpha
#     if(is_tagged){
#       if(n_loc>1){
#         if (!is.na(r_l[i])) { #recapture
#           idx_det <- if (n_loc > 1) n_loc + r_l[i] else n_loc + 1
#           if (idx_det > length(alpha) || idx_det < 1)
#             stop("idx_det out of bounds: ", idx_det)
#           prob_rec <- alpha[idx_det]
#         } else { #tagged but never recovered
#           prob_rec <- sum(alpha[ (n_loc + 1):(2 * n_loc) ])
#         }
#       }else{
#         # if (!is.na(r_l[i])) { #recapture
#           idx_det <- if (n_loc > 1) n_loc + r_l[i] else n_loc + 1
#           if (idx_det > length(alpha) || idx_det < 1)
#             stop("idx_det out of bounds: ", idx_det)
#           prob_rec <- alpha[2]
#         # } else { #tagged but never recovered
#           # prob_rec <- sum(alpha[ (n_loc + 1):(2 * n_loc) ])
#         # }
#       }
#     }
#
#     if (entry == s) prob_rec <- 0         # can’t be seen again after last primary
#
#
#     if (is_tagged) {
#       # These fish were all tagged at first detection
#       nll <- nll - dbinom(count, size = count, prob = p_marked, log = TRUE)
#
#       if (!is.na(recapture)) {
#         nll <- nll - dbinom(count, size = count, prob = prob_rec, log = TRUE)
#       } else {
#         nll <- nll - dbinom(0, size = count, prob = prob_rec, log = TRUE)
#       }
#
#     } else {
#       # These fish were detected but not tagged — they all failed to be tagged
#       nll <- nll - dbinom(count, size = count, prob = 1 - p_marked, log = TRUE)
#     }
#   }
#
#   if(length(struct$Z_list$phi)>0){
#     nll <- nll - sum(dnorm(u_phi,0,exp(log_sd_phi),TRUE))
#     out$u_phi <- u_phi
#   }
#   if(length(struct$Z_list$p)>0){
#     nll <- nll - sum(dnorm(u_p,0,exp(log_sd_p),TRUE))
#     out$u_p <- u_p
#   }
#
#   out$p_for_phi_group <- p_for_phi_group  <- map_groups(phi_group, p_group)
#   out$w_for_phi_group <- w_for_phi_group  <- map_groups(phi_group, w_group)
#
#   calc_psi <- function(S, D, pent_mat, sp_pr, s, n_loc)
#   {
#
#     G   <- dim(S)[1]
#     n_states  <- 2 * n_loc + 1
#     psiPtot_g <- numeric(G)
#     detect_mat<- matrix(0, G, s)
#
#     if (n_loc == 1)
#       sp_pr <- 1
#     else
#       sp_pr <- as.vector(sp_count$n / sum(sp_count$n))
#
#     for (g in 1:G) {
#
#       alpha <- numeric(n_states)
#       alpha[1:n_loc] <- pent_mat[g,1] * sp_pr          # week‑1 entrants
#
#       # gp <- p_group_for_phi_group[g]
#       gp <- p_for_phi_group[g]
#       # gw <- w_for_phi_group[g]
#       # print(c(gp,gw))
#       for (t in 1:s) {
#
#         ## continuous‑time entrants during week t
#         for (j in 1:n_loc)
#           alpha[j] <- alpha[j] +
#             pent_mat[1,t] * sp_pr[j] *
#             (S[g, j, j] - 1) / log(S[g, j, j])
#
#         ## ---------- FIX: use off‑diagonal p, not diagonal -------------
#         p_vec       <- D[gp , 1:n_loc , 1:n_loc + n_loc]
#         detect_g_t  <- sum(alpha[1:n_loc] * p_vec)
#         detect_mat[g, t] <- detect_g_t
#
#         ## propagate one interval
#         alpha <- alpha %*% S[g,,] %*% D[gp,,]
#       }
#
#       psiPtot_g[g] <- sum(detect_mat[g, ])
#     }
#
#     list(psiPtot_g = psiPtot_g,
#          detect_mat = detect_mat)
#   }
#   #
#   out_g <- calc_psi(
#     S      = S,          # G × (2n+1) × (2n+1)
#     D  = D,      # G × n_loc   (or 1 × n_loc if pooled)
#     pent_mat   = pent_mat,
#     sp_pr  = sp_pr,
#     s      = s,
#     n_loc  = n_loc)
#   #
#   psiPtot_g  <- out_g$psiPtot_g    # length G
#   detect_mat <- out_g$detect_mat   # G × s
#   out$psiPtot_g <- out_g$psiPtot_g
#   out$detect_mat <- out_g$detect_mat
#
#   # survival (phi / p) varies by YEAR  ->  G levels
#   ## ---------- helper -----------------------------------------------
#   # get_phi_groups <- function(struct)
#   # {
#   #   col <- struct$groups$phi$group_phi_re1
#   #   if (is.null(col) || length(col) == 0)                  # <- no groups
#   #     factor(rep(1L, length(data$t_k)))                    # one dummy level
#   #   else
#   #     factor(col)
#   # }
#
#   if(!is.null(struct$groups$phi$group_phi_re1)){
#     year_fac <- G <- max(struct$groups$phi$group_phi_re1)         # now == 1
#     out$G <- G
#   }else{
#     year_fac <- G <- 1         # now == 1
#     out$G <- G
#   }
#
#   ## K is handled the same way for Nsuper
#   get_ns_groups <- function(struct) {
#     col <- struct$groups$Nsuper[[1]]
#     if (is.null(col) || length(col) == 0) factor(1L) else factor(col)
#   }
#   ns_fac <- get_ns_groups(struct)
#   K      <- nlevels(ns_fac)
#
#   # abundance N* is to be estimated by NS_GROUP  ->  K levels
#   if(length(struct$groups$Nsuper[[1]])==0){
#     ns_fac <- factor(1L)
#     K <- nlevels(ns_fac)
#     out$K <- K
#   }else{
#     ns_fac       <- factor(struct$groups$Nsuper[[1]])        # same length as year_fac
#     K            <- nlevels(ns_fac)
#     out$K <- K
#   }
#
#   #
#   A <- matrix(0L, nrow = K, ncol = G)
#   A[cbind(as.integer(ns_fac), as.integer(year_fac))] <- 1L
#
#
#   ## for every survival‑group g, find its Ns_group k and put a 1
#   for (g in seq_len(G)) {
#     year_label <- levels(year_fac)[g]          # e.g. "2018"
#     k          <- which(levels(ns_fac) ==
#                           ns_fac[ match(year_label, year_fac) ])
#     A[k, g] <- 1L
#   }
#   out$A <- A
#
#   # option 1 – convert the numeric matrix to an AD matrix once
#   A_ad      <- RTMB::matrix(A)          # keeps dimensions and becomes AD type
#
#   print('A')
#   print(A)
#   print('psiPtot_g')
#   print(psiPtot_g)
#   ## make psiPtot_g an AD-matrix with K rows, 1 column ------------
#   psiPtot_mat <- RTMB::matrix(psiPtot_g, nrow = length(psiPtot_g), ncol = 1)
#   ## matrix product:  (K × K)  %*%  (K × 1)  →  (K × 1)
#
#   psi_dag <- A %*% psiPtot_g
#   out$psi_dag <- psi_dag
#
#   print('psi_dag')
#   print(psi_dag)
#   if (length(struct$Z_list$Nsuper) > 0L) {          # with RE
#     ## fixed part recycled to length K, then add the random effect
#     eta_N  <- rep(beta_Nsuper[1L], K) + u_Nsuper     # linear predictor (log-scale)
#
#     Nsuper_vec <- exp(eta_N)                        # length K
#     lambda_k   <- Nsuper_vec * psi_dag              # mean of Poisson
#     ## prior for the random effects
#     nll <- nll - sum(dnorm(u_Nsuper, 0, exp(log_sd_Nsuper), log = TRUE))
#   } else {                                          # no RE
#     Nsuper_vec <- rep(exp(beta_Nsuper[1L]), K)      # same scalar for every k
#     lambda_k   <- Nsuper_vec * psi_dag
#   }
#   out$Nsuper_vec <- Nsuper_vec
#
#   nll <- nll - sum(dpois(uTot_g, lambda_k, log = TRUE))
#   out$uTot_g <- uTot_g
#   out$lambda_k <- lambda_k
#
#   u_k_t <- matrix(0, K, s)
#   for (i in seq_along(t_k))
#     u_k_t[ struct$groups$Nsuper[[1]][i] , d$t_k[i] ] <-
#     u_k_t[ struct$groups$Nsuper[[1]][i] , d$t_k[i] ] + d$n[i]
#
#   out$detect_k_t <- detect_k_t <- A %*% detect_mat          # K × s
#   out$multP_k_t <- multP_k_t  <- detect_k_t / rowSums(detect_k_t)
#
#   for (k in 1:K)
#     nll <- nll - dmultinom(u_k_t[k, ], prob = multP_k_t[k, ], log = TRUE)
#
#
#   RTMB::REPORT(out)
#   return(nll)
#
# }
