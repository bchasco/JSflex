rm(list=ls())
set.seed(1)
df <- data.frame(
  f_tk   = factor(rep(1:5, each = 30)),
  f_tl   = factor(sample(c("A", "B", "C"), 50, replace = TRUE)),
  f_yr   = factor(sample(2010:2012, 50, replace = TRUE)),
  f_site = factor(sample(1:3, 50, replace = TRUE))
)

formulas <- list(
  phi = ~ -1 + f_yr + f_tk + (1|f_site),
  p   = ~ 1
)

design <- make_design_list_2(formulas, list(state = NULL, input = "f_tk"), df)
lk <- make_process_lookup(c("phi"), design, df, state_var = design$state, period_var = "f_tk")
head(lk$phi$lookup)
table(lk$phi$lookup$state)

params    <- make_param_vectors(design, c("phi"))

# fixed effects get 1,2,3…
beta_vals <- seq_along(params$beta$phi)
names(beta_vals) <- names(params$beta$phi)
# random effects get 100,200,300…
u_vals    <- seq_along(params$u$phi) * 100
names(u_vals)   <- names(params$u)


# (1) precompute the LP for each minimal‐design row
eta_min <- as.numeric(
  (lk$phi$X_min) %*% beta_vals +
     lk$phi$Z_min %*% u_vals
)
eta_row <- as.numeric(
  rownames(lk$phi$X_min) #%*% beta_vals +
  # design$Z_group$phi #%*% u_vals
)

# dims
G <- max(lk$phi$lookup$g)
S <- max(lk$phi$lookup$state) + 2  # +2 for the “dead” and “unseen” states
T <- max(lk$phi$lookup$time)
S2 <- array(0, c(G,S,S,T))
Sidx <- array(0, c(G,S,S,T))

# eta_min <- result$phi$X_min %*% beta_vals    # + result$phi$Z_min %*% u_vals
idx     <- cbind(lk$phi$lookup$g, lk$phi$lookup$state, lk$phi$lookup$state, lk$phi$lookup$time)
S2[idx] <- eta_min[ lk$phi$lookup$design_row ]
Sidx[idx] <- eta_row[ lk$phi$lookup$design_row ]
S2[1,1,,]
Sidx[1,1,,]
S2[1,1,,]
S2[1,S-2,,]
S2[G,S-2,,]
print(dim(S2))
