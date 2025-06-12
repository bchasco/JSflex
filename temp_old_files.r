#' #' Build a grouping factor by dropping any 'period' term (e.g., f_tk)
#' #'
#' #' @param fml A one- or two-sided formula with fixed and/or random effects.
#' #' @param data A data.frame with all variables in the formula.
#' #' @param period Name of the time/period variable to exclude (default = "f_tk").
#' #'
#' #' @return A factor representing interaction of all non-period grouping variables.
#' make_group_without_period <- function(fml, data, period = "f_tk", state_var) {
#'   # Extract all RHS variables (fixed + random)
#'   bars <- lme4::findbars(fml)
#'   re_vars <- if (length(bars)) all.vars(bars[[1]][[3]]) else character(0)
#'
#'   rhs <- if (length(fml) == 2) fml[[2]] else fml[[3]]
#'   fe_vars <- all.vars(rhs)
#'
#'   # Combine and drop period variable
#'   if(!is.null(state_var)){
#'     exclude_vars <- c('f_tl',period)
#'     print("test")
#'   }else{
#'     exclude_vars <- c(period)
#'   }
#'   vars <- setdiff(unique(c(re_vars, fe_vars)), exclude_vars)
#'
#'   # Return interaction of all remaining variables (or single group if none)
#'   if (length(vars) == 0L) {
#'     factor(rep(1L, nrow(data)))
#'   } else {
#'     interaction(data[vars], drop = TRUE)
#'   }
#' }
#'
#' #' Default grouping helper
#' #'
#' #' Combine any factor columns in a model frame and random-effect factors
#' #' to produce a single grouping factor.
#' #'
#' #' @param mf_full A model.frame used to extract fixed-effect factors.
#' #' @param flist A list of random-effect factor vectors.
#' #'
#' #' @return A factor representing the interaction of all group factors.
#' make_process_group <- function(mf_full, flist) {
#'   fixed_fac <- Filter(is.factor, mf_full)
#'   rand_fac  <- if (length(flist)) as.data.frame(flist) else mf_full[ , 0, drop = FALSE]
#'   if (ncol(fixed_fac) + ncol(rand_fac) == 0) {
#'     factor(rep(1L, nrow(mf_full)))
#'   } else {
#'     interaction(cbind(fixed_fac, rand_fac), drop = TRUE)
#'   }
#' }
#'
#' #' Extract design matrices and grouping structures
#' #'
#' #' Build fixed-effect (X) and random-effect (Z) design matrices and group factors
#' #' for a list of model components.
#' #'
#' #' @param formulas A named list of formulas (one per process).
#' #' @param data A data.frame with covariates.
#' #'
#' #' @return A list with X_list, Z_list, Z_terms, groups, re_factor, and map_to_joint.
#' make_design_list <- function(formulas, state, data) {
#'   data <- as.data.frame(data)
#'   data$._dummy_response <- 0L
#'   design_list <- list(
#'     X_list       = list(),
#'     Z_list       = list(),
#'     Z_terms      = list(),
#'     groups       = list(),
#'     group_factor = list(),
#'     re_factor    = list(),
#'     group_noperiod = list()
#'   )
#'   for (proc in names(formulas)) {
#'     fml <- formulas[[proc]]
#'     if (length(fml) == 2)
#'       fml <- reformulate(deparse(fml[[2]]), response = "._dummy_response")
#'     mf_full <- model.frame(lme4::subbars(fml), data)
#'     design_list$X_list[[proc]] <- model.matrix(lme4::nobars(fml), mf_full)
#'     bars <- lme4::findbars(fml)
#'     Zt_trm <- if (length(bars)) lme4::mkReTrms(bars, mf_full) else NULL
#'     if (!is.null(Zt_trm)) {
#'       design_list$Z_terms[[proc]] <- Zt_trm
#'       design_list$Z_list[[proc]]  <- t(as.matrix(Zt_trm$Zt))
#'       design_list$groups[[proc]]   <-
#'         setNames(lapply(Zt_trm$flist, as.integer),
#'                  paste0("group_", proc, "_re", seq_along(Zt_trm$flist)))
#'       design_list$re_factor[[proc]] <- factor(rep(1L, nrow(mf_full)))
#'     } else {
#'       design_list$Z_terms[[proc]] <- NULL
#'       design_list$Z_list[[proc]]  <- matrix(0L, nrow(mf_full), 0L)
#'       design_list$groups[[proc]]   <- list()
#'       design_list$re_factor[[proc]] <- factor(rep(1L, nrow(mf_full)))
#'     }
#'     design_list$group_factor[[proc]] <-
#'       make_process_group(mf_full, if (is.null(Zt_trm)) list() else Zt_trm$flist)
#'   }
#'   for (proc in intersect(names(formulas), c("w","phi","p","t_var","Nsuper"))) {
#'     design_list$group_noperiod[[proc]] <-
#'       make_group_without_period(formulas[[proc]], data, state_var = state)
#'   }
#'   design_list
#' }
#'
#' #' Map joint factor levels to a process factor
#' #'
#' #' @param joint_fac Factor of joint grouping.
#' #' @param proc_fac Factor of process grouping.
#' #'
#' #' @return Integer vector mapping each joint level to a proc level.
#' map_joint_to_proc <- function(joint_fac, proc_fac) {
#'   sapply(levels(joint_fac), function(lvl) {
#'     as.integer(proc_fac[ match(TRUE, joint_fac == lvl) ])
#'   })
#' }
#'
#' #' Initialize parameters from a design list
#' #'
#' #' @param design Output of make_design_list().
#' #' @param init_sd Initial standard deviation for random effects.
#' #'
#' #' @return Named list of starting parameters (beta_*, u_*, log_sd_*).
#' make_all_params_2.1 <- function(design, data, init_sd = 0.1) {
#'   params <- list()
#'   for (proc in names(design$group_factor)) {
#'     X <- design$X_list[[proc]]; Z <- design$Z_list[[proc]]
#'     params[[paste0("beta_", proc)]] <- if (ncol(X) > 0) rep(0, ncol(X)) else NA
#'     if (ncol(Z) > 0) {
#'       params[[paste0("u_", proc)]]      <- rep(0, ncol(Z))
#'       params[[paste0("log_sd_", proc)]] <- log(init_sd)
#'     } else {
#'       params[[paste0("u_", proc)]]      <- NA
#'       params[[paste0("log_sd_", proc)]] <- NA
#'     }
#'     if(!is.null(data$state) & proc == "t_var"){
#'       params[[paste0("u_", proc)]]      <- 0
#'       params[[paste0("log_sd_", proc)]] <- 1
#'     }
#'   }
#'
#'   params$beta_Nsuper <- rep(10, length(params$beta_Nsuper))
#'   params$beta_phi <- rep(-2,length(params$beta_phi))
#'   params$beta_p <- rep(-1,length(params$beta_p))
#'
#'   if(!is.null(data$state)){
#'     n_states <- 5
#'     n_par_M <- n_states * (n_states - 1) / 2
#'     params$u_t_var <- rep(0,n_par_M)
#'   }
#'
#'
#'   params
#' }
#'
#' #' Build RTMB data list from design object and data
#' #'
#' #' This function centralizes and simplifies the construction of the full RTMB data list.
#' #' It pulls together the core fields, grouping structures, joint mappings, and scalar values.
#' #'
#' #' @param design A design list from `make_design_list()`
#' #' @param data   A data.frame with individual-level capture data
#' #' @param state  Optional. Name of column indicating state (e.g., `"t_l"`)
#' #' @param period Character. Column name for the period variable (default: "t_k")
#' #'
#' #' @return A named list suitable for RTMB model fitting.
#' #' @export
#' make_RTMB_data_list <- function(design, data, state = NULL, period = "t_k") {
#'
#'   phi_group <- design$group_noperiod[['phi']]
#'   p_group   <- design$group_noperiod[['p']]
#'   w_group   <- design$group_noperiod[['w']]
#'   tvar_group <- design$group_noperiod[['t_var']]
#'   Nsuper_group <- design$group_noperiod[['Nsuper']]
#'
#'   joint_group <- interaction(phi_group, p_group, w_group, tvar_group, drop = TRUE)
#'   G <- nlevels(joint_group)
#'
#'   phi_for_joint <- map_joint_to_proc(joint_group, phi_group)
#'   p_for_joint   <- map_joint_to_proc(joint_group, p_group)
#'   w_for_joint   <- map_joint_to_proc(joint_group, w_group)
#'   tvar_for_joint <- map_joint_to_proc(joint_group, tvar_group)
#'   Nsuper_for_joint <- map_joint_to_proc(joint_group, Nsuper_group)
#'
#'   if (nlevels(Nsuper_group) == 1) {
#'     K <- 1L
#'     A <- matrix(1, nrow = 1, ncol = G)
#'   } else {
#'     K <- nlevels(Nsuper_group)
#'     A <- Matrix::sparseMatrix(
#'       i = as.integer(Nsuper_group),
#'       j = as.integer(joint_group),
#'       x = 1,
#'       dims = c(K, G)
#'     )
#'     A <- as.matrix(A)
#'   }
#'
#'   # Construct u (number tagged per time step)
#'   s <- max(na.omit(c(data$t_k, data$r_k)))
#'   u_vec <- numeric(s)
#'   u_sum <- aggregate(list(n = data$n), by = list(t = data[[period]]), sum)
#'   u_vec[u_sum$t] <- u_sum$n
#'
#'   # Spatial proportions of detections
#'   sp_count <- aggregate(list(n = data$n), by = list(loc = data[[state %||% 't_l']]), sum)
#'   sp_pr <- sp_count$n / sum(sp_count$n)
#'
#'   # uTot_g per Nsuper_group
#'   if (length(Nsuper_group) > 1) {
#'     uTot_g <- as.numeric(tapply(data$n, Nsuper_group, sum))
#'   } else {
#'     uTot_g <- sum(data$n)
#'   }
#'
#'   list(
#'     t_k     = data$t_k,
#'     r_k     = data$r_k,
#'     t_l     = data$t_l,
#'     f_tl     = data$f_tl,
#'     r_l     = data$r_l,
#'     tag     = as.integer(data$tag),
#'     n       = data$n,
#'     s       = s,
#'     period_fac = data[[period]],
#'     X_list  = design$X_list,
#'     Z_list  = design$Z_list,
#'     group_factor = design$group_factor,
#'     phi_for_joint = phi_for_joint,
#'     p_for_joint   = p_for_joint,
#'     w_for_joint   = w_for_joint,
#'     tvar_for_joint = tvar_for_joint,
#'     Nsuper_for_joint = Nsuper_for_joint,
#'     joint_group = joint_group,
#'     Nsuper_group = Nsuper_group,
#'     G = G,
#'     K = K,
#'     A = A,
#'     u = u_vec,
#'     uTot_g = uTot_g,
#'     sp_pr = sp_pr,
#'     state = state
#'   )
#' }
#'
#' make_group_mapping_df <- function(joint_group, group_factor) {
#'   G <- nlevels(joint_group)
#'   data.frame(
#'     joint = levels(joint_group),
#'     lapply(group_factor, function(gf) {
#'       vapply(levels(joint_group), function(lvl) {
#'         idx <- which(joint_group == lvl)[1]
#'         as.integer(gf[idx])
#'       }, integer(1))
#'     })
#'   )
#' }
#'
#'
plot_Sankey <- function(data) {
  library(viridis)
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(scales)

  mapping_df <- data.frame(
    joint = levels(data$joint_group),
    phi   = data$phi_for_joint,
    p     = data$p_for_joint,
    w     = data$w_for_joint,
    t_var = data$tvar_for_joint,
    Nsuper = data$Nsuper_for_joint
  ) %>%
    pivot_longer(cols = -joint, names_to = "process", values_to = "group") %>%
    mutate(
      joint   = factor(joint, levels = unique(joint)),
      process = factor(process, levels = c("phi", "p", "w", "t_var","Nsuper")),
      group   = as.integer(group)
    ) %>%
    group_by(process) %>%
    mutate(group_scaled = rescale(group, to = c(1, 100))) %>%  # Normalize
    ungroup()

  ggplot(mapping_df, aes(x = process, y = joint)) +
    geom_tile(aes(fill = group_scaled), color = "white") +
    geom_text(aes(label = group), size = 2.5, color = "black") +
    scale_fill_viridis_c(name = "Group", option = "plasma") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 7),
      panel.grid = element_blank()
    ) +
    labs(
      title = "Group Mapping by Process (Normalized Color by Process)",
      x = "Process",
      y = "Joint Group"
    )
  # invisible(mapping_df)
}

# make_random <- function(params) {
#   re_names <- c("u_phi", "u_p", "u_w", "u_v", "u_Nsuper", "u_t_var")
#   random <- names(Filter(function(x) length(x) > 1, params[re_names]))
#   return(random)
# }
#'
#' make_map <- function(params) {
#'   map <- list()
#'
#'   # Always fixed (you can customize this logic)
#'   map$beta_v     <- factor(NA)
#'   map$beta_w     <- factor(NA)
#'   map$u_v        <- factor(NA)
#'   map$log_sd_v   <- factor(NA)
#'   map$u_Nsuper        <- factor(NA)
#'   map$log_sd_Nsuper   <- factor(NA)
#'   map$beta_t_var <- factor(NA)
#'   if(is.na(params$log_sd_t_var)){
#'     map$u_t_var <- factor(rep(NA,length(params$u_t_var)))  # optionally comment/uncomment
#'     map$log_sd_t_var <- factor(NA)
#'   }
#'   if(is.na(params$log_sd_phi)){
#'     map$u_phi <- factor(rep(NA,length(params$u_phi)))  # optionally comment/uncomment
#'     map$log_sd_phi <- factor(NA)
#'   }
#'   if(is.na(params$log_sd_p)){
#'     map$u_p <- factor(rep(NA,length(params$u_p)))  # optionally comment/uncomment
#'     map$log_sd_p <- factor(NA)
#'   }
#'   if(is.na(params$log_sd_w)){
#'     map$u_w <- factor(rep(NA,length(params$u_w)))  # optionally comment/uncomment
#'     map$log_sd_w <- factor(NA)
#'   }
#'
#'
#'   return(map)
#' }
