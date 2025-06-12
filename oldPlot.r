## helper to stack a matrix to long format  ------------------------
mat_long <- function(mat, var = "value") {
  as.data.frame(mat) |>
    tibble::rownames_to_column("grp") |>
    tidyr::pivot_longer(-grp, names_to = "t", values_to = var) |>
    dplyr::mutate(t   = as.numeric(gsub("^V", "", t)),
                  grp = as.factor(grp))
}

## ----------------------------------------------------------------
## 1.  first‑detection fit: observed counts vs. model mean ---------
## ----------------------------------------------------------------
obs <- data.frame(
  t   = data$t_k,
  n   = data$n,
  grp = factor(d$f_sex, levels = levels(d$f_sex))       # or Nsuper group
) |>
  dplyr::group_by(t, grp) |>
  dplyr::summarise(n = sum(n), .groups = "drop")

## ---- expected counts = per‑fish p  ×  N*  -----------------------
Nvec   <- rep$out$Nsuper_vec                       # length K (or G)
lambda <- sweep(rep$out$detect_k_t, 1, Nvec, "*")   # rows × cols

pred <- mat_long(lambda, "exp_n")
levels(pred$grp) <- levels(d$f_sex)

## ---- plot -------------------------------------------------------
library(ggplot2)
ggplot(obs, aes(t, n)) +
  geom_col(fill = "grey80") +
  geom_line(data = pred, aes(t, exp_n), colour = "red", linewidth = 1.1) +
  facet_wrap(~ grp) +
  labs(title = "Observed first detections vs. model expectation",
       y = "count", x = "survey occasion") +
  theme_bw()

print(rep$out$Nsuper_vec)

# #
# ## -----------------------------------------------------------------
# ## 4. quick residual plot for the CJS part --------------------------
# ## -----------------------------------------------------------------
# cjs_df <- data.frame(
#   p_marked = rep$out$v_w[data$t_k],
#   n        = data$n
# ) |>
#   dplyr::mutate(
#     exp = n * p_marked,
#     res = n - exp
#   )
#
# ggplot(cjs_df, aes(exp, res)) +
#   geom_point(alpha = 0.4) +
#   geom_hline(yintercept = 0, linetype = 2, colour = "red") +
#   labs(title = "Residuals for first‑detection (CJS block)",
#        x = "expected", y = "observed – expected")
# ## -----------------------------------------------------------------
# ## 2. entry‐timing distribution (pent) ------------------------------
# ## -----------------------------------------------------------------
# pent_df <- data.frame(
#   t   = seq_along(rep$out$pent_mat),
#   pent_mat = t(rep$out$pent_mat)
# )
#
# ggplot(pent_df, aes(t, pent_mat)) +
#   geom_col(fill = "steelblue") +
#   labs(title = "Estimated temporal entry distribution (pent)",
#        x = "survey occasion", y = "probability")
#
# #
# ## ------------------------------------------------------
# ## 1.  SURVIVAL ϕ  (one value per φ‑group)  -------------
# ## ------------------------------------------------------
# n_loc <- 1                               # change if >1
#
# phi_g <- rep$out$S[, 1, 1]                   # S[g, live, live]
# df_phi <- data.frame(
#   grp = seq_along(phi_g),
#   phi = phi_g
# )
#
# ggplot(df_phi, aes(factor(grp), phi)) +
#   geom_point(size = 3) +
#   labs(title = "Daily survival probability (ϕ)",
#        x = "ϕ‑group", y = "ϕ") +
#   ylim(0, 1)
# #
# #
# # ## ------------------------------------------------------
# # ## 2.  DETECTION p  (one value per p‑group) -------------
# # ## ------------------------------------------------------
# # p_g <- fit$rep$out$D[, 1, 1 + n_loc]             # D[g, live, detected]
# # df_p <- data.frame(
# #   grp = seq_along(p_g),
# #   p   = p_g
# # )
# #
# # ggplot(df_p, aes(factor(grp), p)) +
# #   geom_point(size = 3, colour = "darkred") +
# #   labs(title = "Detection probability (p)",
# #        x = "p‑group", y = "p") +
# #   ylim(0, 1)
# #
# #
# # ## ------------------------------------------------------
# # ## 3.  N* (super‑population size)  ----------------------
# # ## ------------------------------------------------------
# # N_k <- fit$rep$out$Nsuper                        # length‑K vector
# # df_N <- data.frame(
# #   grp = seq_along(N_k),
# #   N   = N_k
# # )
# #
# # ggplot(df_N, aes(factor(grp), N)) +
# #   geom_col(fill = "grey70") +
# #   labs(title = "Estimated super‑population size (N*)",
# #        x = "N*‑group", y = "individuals")
#
# print(rep$out$Nsuper_vec)
