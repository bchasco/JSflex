rm(list=ls())
set.seed(1)
df <- data.frame(
  f_tk   = factor(rep(1:5, each = 30)),
  f_tl   = factor(sample(c("A", "B", "C"), 50, replace = TRUE)),
  f_yr   = factor(sample(2010:2012, 50, replace = TRUE)),
  f_site = factor(sample(1:3, 50, replace = TRUE))
)

formulas <- list(
  phi = ~ -1 + (1|f_tl) + f_tk,
  p   = ~ 1
)

design <- make_design_list_2(formulas, list(state = 'f_tl', input = "f_tk"), df)
lk <- make_process_lookup(c("phi"), design, df, state_var = design$state, period_var = "f_tk")
params    <- make_param_vectors(design, c("phi"))

# fixed effects get 1,2,3…
beta_vals <- seq_along(params$beta$phi)
names(beta_vals) <- names(params$beta$phi)
# random effects get 100,200,300…
if(ncol(lk$phi$Z_min)>0){
  u_vals    <- seq_along(params$u$phi) * 100
} else {
  u_vals <- NA
}

G <- max(lk$phi$lookup$g)
S <- if (!is.null(design$state)) nlevels(df$f_tl) + 2 else 3L
T <- max(lk$phi$lookup$time)

S2 <- array(0, c(G, S, S, T))


# for(i in 1:nrow(lk$phi$lookup)){
#   print(i)
#   row <- lk$phi$lookup[i,]
#   beta_i <- lk$phi$X_min[row$design_row,] %*% beta_vals
#   if(ncol(lk$phi$Z_min)>0){
#     re_i <- lk$phi$Z_min[row$re_row,] %*% u_vals
#     beta_i <- beta_i + re_i
#   } else {
#     re_i <- 0
#   }
#   S2[row$g,row$state,row$state,row$time] <- beta_i
# }
# S2[1,1,,1]
# S2[1,,,1]
#
# S2[1,1,,]
# S2[1,S-2,,]
# S2[G,S-2,,]
# S2[1,,,1]
print(dim(S2))
table(S2)
