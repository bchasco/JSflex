#' @title Build Lookup Tables from Group-State-Time to Design Row for Multiple Processes
#' @description
#' For one or more model processes, generates lookup tables mapping combinations of
#' group, state, and time to row indices in the minimal design matrices. Necessary
#' for populating full S2 or D2 arrays even when some combinations are unobserved.
#'
#' @param procs Character vector of process names (e.g., c("phi","p")).
#' @param design Output of \code{make_design_list_2()}.
#' @param data The original dataset.
#' @param state_var Optional string name of the state variable (e.g., "f_tl"); if \code{NULL}, a single alive state is assumed.
#' @param period_var String name of the time (period) variable (e.g., "f_tk").
#'
#' @return A named list of length \code{length(procs)}, each element a list with:
#' \describe{
#'   \item{lookup}{data.frame with columns \code{g}, \code{state}, \code{time}, \code{design_row}, \code{re_row}}
#'   \item{X_min}{Minimal fixed-effect design matrix for that process}
#'   \item{Z_min}{Minimal random-effect design matrix for that process (if applicable)}
#' }
#'
#' @examples
#' df <- data.frame(
#'   f_tk   = rep(1:5, each = 6),
#'   f_tl   = factor(sample(c("A", "B", "C"), 30, TRUE)),
#'   f_yr   = factor(sample(2010:2012, 30, TRUE)),
#'   f_site = factor(sample(1:3, 30, TRUE))
#' )
#'
#' # 1a) fixed period effect
#' formulas1 <- list(phi = ~ -1 + f_yr + f_tk)
#' design1 <- make_design_list_2(formulas1, list(state="f_tl", input="f_tk"), df)
#' res1 <- make_process_lookup("phi", design1, df, state_var="f_tl", period_var="f_tk")
#'
#' # 1b) random period effect
#' formulas2 <- list(phi = ~ -1 + f_yr + (1|f_tk))
#' design2 <- make_design_list_2(formulas2, list(state="f_tl", input="f_tk"), df)
#' res2 <- make_process_lookup("phi", design2, df, state_var="f_tl", period_var="f_tk")
#'
#' # 1c) interaction period effect
#' formulas3 <- list(phi = ~ -1 + f_yr:f_tk)
#' design3 <- make_design_list_2(formulas3, list(state="f_tl", input="f_tk"), df)
#' res3 <- make_process_lookup("phi", design3, df, state_var="f_tl", period_var="f_tk")
#'
#' # multiple processes at once
#' formulas_all <- list(phi=~f_yr+f_tl, p=~1)
#' design_all <- make_design_list_2(formulas_all, list(state="f_tl", input="f_tk"), df)
#' res_all <- make_process_lookup(c("phi","p"), design_all, df, state_var="f_tl", period_var="f_tk")
#' str(res_all)
#'
#' @export
make_process_lookup <- function(procs, design, data,
                                state_var  = NULL,
                                period_var = NULL) {
  procs <- as.character(procs)
  out   <- setNames(vector("list", length(procs)), procs)

  for (proc in procs) {
    fml <- design$formulas[[proc]]
    txt <- paste(deparse(fml), collapse=" ")
    has_time  <- !is.null(period_var) && grepl(paste0("\\b",period_var,"\\b"), txt)
    has_state <- !is.null(state_var)  && grepl(paste0("\\b", state_var,"\\b"), txt)

    # 1) choose which grouping defines G
    g_fac <- switch(
      paste0(has_time,"_",has_state),
      "FALSE_FALSE" = design$group_excl_time_state[[proc]],
      "TRUE_FALSE"  = design$group_noperiod        [[proc]],
      "FALSE_TRUE"  = design$group_excl_time_state [[proc]],
      "TRUE_TRUE"   = design$group_factor          [[proc]]
    )
    n_g <- nlevels(g_fac)

    # 2) pick exactly one row per G
    rep_idx <- match(levels(g_fac), g_fac)
    dd0     <- data[rep_idx, , drop=FALSE]

    # 3) build the full G × state × time grid
    n_state <- if (has_state) nlevels(data[[state_var]]) else 1L
    time_vals <- if (has_time) sort(unique(data[[period_var]])) else 1
    n_time <- length(time_vals)

    grid <- expand.grid(
      g     = seq_len(n_g),
      state = seq_len(n_state),
      time  = seq_len(n_time)
    )

    # 4) expand dd0 to match that grid, *overwriting* state/period
    dd <- dd0[ grid$g, , drop=FALSE ]
    if (has_state) dd[[state_var]]  <- levels(data[[state_var]])[ grid$state ]
    if (has_time ) dd[[period_var]] <-   time_vals[      grid$time  ]

    # 5) rebuild minimal design on that expanded data
    mf2   <- model.frame(lme4::subbars(fml), dd)
    X_min <- model.matrix( lme4::nobars(fml), mf2 )

    bars <- lme4::findbars(fml)
    if (length(bars)) {
      Zt2  <- lme4::mkReTrms(bars, mf2)
      Z_min <- t(as.matrix(Zt2$Zt))
    } else {
      Z_min <- matrix(0, nrow(X_min), 0)
    }

    # 6) now you have perfect 1–1 mapping
    grid$design_row <- seq_len(nrow(grid))
    grid$re_row     <- if (ncol(Z_min) > 0) grid$design_row else NA_integer_

    out[[proc]] <- list(
      lookup = grid,
      X_min  = X_min,
      Z_min  = Z_min
    )
  }

  out
}

#' @title Create Parameter Vectors for Process
#' @description
#' For a given process name, extract the number of fixed and random effect parameters
#' based on the design object and create named beta and u parameter vectors.
#'
#' @param design Output of `make_design_list_2()`.
#' @param proc A string giving the process name (e.g., "phi" or "p").
#'
#' @return A list with two elements:
#' - `beta`: Named vector of fixed effects for this process.
#' - `u`: Named vector of random effects for this process (may be empty).
#'
#' @examples
#' df <- data.frame(
#'   f_tk   = rep(1:5, each = 6),
#'   f_tl   = factor(sample(c("A", "B", "C"), 30, replace = TRUE)),
#'   f_yr   = factor(sample(2010:2012, 30, replace = TRUE)),
#'   f_site = factor(sample(1:3, 30, replace = TRUE))
#' )
#'
#' formulas <- list(
#'   phi = ~ f_yr + (1|f_site),
#'   p   = ~ 1
#' )
#'
#' design <- make_design_list_2(formulas, list(state = "f_tl", input = "f_tk"), df)
#' make_param_vectors(design, "phi")
#' make_param_vectors(design, "p")
#'
#' @export
make_param_vectors <- function(design, processes = names(design$formulas)) {
  beta <- list()
  u    <- list()

  for (proc in processes) {
    X <- design$X_group[[proc]]
    beta[[proc]] <- setNames(rep(0, ncol(X)), colnames(X))

    Z <- design$Z_group[[proc]]
    if (!is.null(Z) && ncol(Z) > 0) {
      u[[proc]] <- setNames(rep(0, ncol(Z)), paste0("u", seq_len(ncol(Z))))
    } else {
      u[[proc]] <- numeric(0)
    }
  }

  list(beta = beta, u = u)
}

#' @title Prepare Group-Based Predictors and Group-to-State Maps for Model Processes
#'
#' @description
#' General utility for extracting linear predictors and group-to-location mappings
#' for any model process (`phi`, `p`, `w`, etc.). It computes:
#' - the linear predictor (`eta`) for each unique group
#' - the mapping from each group to one or more states (locations)
#' - a complete index grid over group × state × period combinations
#'
#' The function ensures that:
#' - the **period variable** (e.g., `f_tk`) never contributes to the grouping structure
#' - the **state variable** (e.g., `f_tl`) is excluded from grouping if provided,
#'   though still used for group-to-location mapping
#'
#' This precomputes the structure required to fully populate matrices like `S2` and `D2`
#' even in the presence of missing group × state × period combinations.
#'
#' @param proc String name of the model process (`"phi"`, `"p"`, etc.)
#' @param data A `data.frame` of covariates, with time and state columns as needed.
#' @param design Output of `make_design_list_2()` (includes formulas and group factors).
#' @param beta_list Named list of fixed effect parameter vectors (`beta$phi`, etc.)
#' @param u_list Named list of random effect parameter vectors (`u$phi`, etc.)
#' @param state_var Optional string naming the state variable (e.g., `"f_tl"`); if `NULL`,
#'   all individuals are assumed to belong to state 1.
#' @param period_var String naming the period variable (default `"f_tk"`).
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{eta}}{Numeric vector of linear predictors for each group.}
#'   \item{\code{g_to_locs}}{List mapping each group to one or more location indices.}
#'   \item{\code{group_fac}}{Factor used to define the unique grouping structure.}
#'   \item{\code{n_groups}}{Number of unique groups.}
#'   \item{\code{mapping}}{A data.frame of all group × state × period combinations.}
#' }
#'
#' @examples
#' # Simulated dataset
#' set.seed(42)
#' df <- data.frame(
#'   f_tk   = rep(1:5, times = 6),
#'   f_tl   = factor(rep(1:2, each = 15)),
#'   f_yr   = factor(rep(2010:2012, length.out = 30)),
#'   f_site = factor(sample(1:3, 30, TRUE))
#' )
#'
#' # Example formulas
#' formulas <- list(
#'   phi = ~ f_yr + (1|f_site),
#'   p   = ~ 1
#' )
#'
#' design <- make_design_list_2(formulas, list(state = "f_tl", input = "f_tk"), data = df)
#' beta  <- list(phi = rep(0.1, ncol(design$X_list$phi)), p = 0)
#' u     <- list(phi = rep(0, ncol(design$Z_list$phi)),   p = numeric(0))
#'
#' param_phi <- prep_param_input("phi", df, design, beta, u, state_var = "f_tl")
#' param_phi$eta         # linear predictor for each group
#' param_phi$g_to_locs   # locations each group covers
#' param_phi$mapping     # all group × state × time combinations
#'
#' @export
prep_param_input <- function(proc, data, design, beta_list, u_list,
                             state_var = NULL, period_var = "f_tk") {
  fml <- design$formulas[[proc]]
  has_period <- grepl(paste0("\\b", period_var, "\\b"),
                      paste(deparse(fml), collapse = " "))

  xz <- get_XZ_group(proc, data, use_noperiod = !has_period)

  grp_fac <- if (!has_period) {
    design$group_noperiod[[proc]]
  } else {
    design$group_factor[[proc]]
  }

  rep_idx <- match(levels(grp_fac), grp_fac)

  loc_idx <- if (!is.null(state_var)) {
    as.integer(data[[state_var]][rep_idx])
  } else {
    rep(1L, length(rep_idx))
  }

  g_to_locs <- split(loc_idx, seq_along(rep_idx))
  eta <- linpred_vec(seq_along(rep_idx), xz, beta_list[[proc]], u_list[[proc]])

  mapping <- expand.grid(
    g = seq_along(rep_idx),
    s = sort(unique(unlist(g_to_locs))),
    t = sort(unique(data[[period_var]]))
  )

  list(
    eta        = eta,
    g_to_locs  = g_to_locs,
    group_fac  = grp_fac,
    n_groups   = length(rep_idx),
    mapping    = mapping
  )
}


#' @title Construct Grouping Factor Excluding a Period Term
#' @description
#' Builds an interaction‐based grouping factor from all fixed‐ and random‐effect
#' variables in a formula, dropping any specified period variable (e.g., `f_tk`)
#' and—if `state_var` is provided—also dropping the state variable `f_tl`.
#'
#' @param fml A one‐ or two‐sided formula with fixed and/or random effects,
#'   e.g. `~ f_yr + (1|f_site)`.
#' @param data A `data.frame` containing all variables in `fml`.
#' @param period A string naming the period variable to exclude (default: `"f_tk"`).
#' @param state_var Optional string naming a state variable; if not `NULL`,
#'   also excludes `"f_tl"`.
#'
#' @return A factor of length `nrow(data)` representing the interaction of all
#'   non‐period (and non‐`f_tl` if `state_var` is used) variables. If no
#'   variables remain, returns a single‐level factor.
#'
#' @examples
#' #--- Setup sample data with both f_tk (time) and f_tl (state) plus other covariates:
#' df <- data.frame(
#'   f_tk  = rep(1:2, each = 6),
#'   f_tl  = rep(c("A","B","C"), times = 4),
#'   f_yr  = factor(rep(2010:2012, length.out = 12)),
#'   f_site = factor(rep(1:3, times = 4))
#' )
#'
#' # 1) Drop only f_tk (period)—group by f_tl, f_yr, f_site:
#' grp_period_only <- make_group_without_period(~ f_tl + f_yr + (1|f_site), df)
#' table(grp_period_only)
#'
#' # 2) Drop only f_tl via state_var—group by f_tk, f_yr, f_site:
#' grp_state_only <- make_group_without_period(~ f_tk + f_yr + (1|f_site),
#'                                            df,
#'                                            period = "f_tk",
#'                                            state_var = "f_tl")
#' table(grp_state_only)
#'
#' # 3) Drop both f_tk and f_tl—group only by f_yr and f_site:
#' grp_both <- make_group_without_period(~ f_yr + (1|f_site),
#'                                      df,
#'                                      period = "f_tk",
#'                                      state_var = "f_tl")
#' table(grp_both)
#'
#' # 4) No other covariates—should return a single level:
#' grp_none <- make_group_without_period(~ (1|f_site),
#'                                      df,
#'                                      period = "f_tk",
#'                                      state_var = "f_tl")
#' table(grp_none)
#'
#' @export
make_group_without_period <- function(fml, data, period = "f_tk", state_var = NULL) {
  bars    <- lme4::findbars(fml)
  re_vars <- if (length(bars)) all.vars(bars[[1]][[3]]) else character(0)

  rhs     <- if (length(fml) == 2) fml[[2]] else fml[[3]]
  fe_vars <- all.vars(rhs)

  exclude_vars <- if (!is.null(state_var)) {
    c("f_tl", period)
  } else {
    period
  }
  vars <- setdiff(unique(c(re_vars, fe_vars)), exclude_vars)

  if (length(vars) == 0L) {
    factor(rep(1L, nrow(data)))
  } else {
    interaction(data[vars], drop = TRUE)
  }
}


#' @title Construct Grouping Factor Excluding Time and (Optionally) State
#' @description
#' Builds an interaction-based grouping factor from all fixed‐ and random‐effect
#' variables in a formula, excluding the time (period) variable (e.g., `f_tk`) and,
#' optionally, a state variable (e.g., `f_tl`). This reduced grouping is essential
#' for determining the number of unique strata in process matrices such as
#' survival, entry, movement, or detection, when time and state are handled separately.
#'
#' @param fml A one‐ or two‐sided formula with fixed and/or random effects,
#'   e.g. `~ f_yr + (1|f_site)`.
#' @param data A `data.frame` containing all variables referenced in `fml`.
#' @param period A string naming the period variable to exclude (default: `"f_tk"`).
#' @param state_var Optional string naming a state variable; if not `NULL`,
#'   also excludes this variable (e.g., `"f_tl"`).
#'
#' @return A factor of length `nrow(data)` representing the interaction of all
#'   fixed and random effect variables in the formula, excluding `period` and
#'   (if specified) `state_var`. If no grouping variables remain, returns a
#'   single-level factor.
#'
#' @examples
#' # Sample data with time (f_tk), state (f_tl), and other covariates
#' df <- data.frame(
#'   f_tk   = rep(1:2, each = 6),
#'   f_tl   = rep(c("A", "B", "C"), times = 4),
#'   f_yr   = factor(rep(2010:2012, length.out = 12)),
#'   f_site = factor(rep(1:3, times = 4))
#' )
#'
#' # 1) Drop only f_tk — group by f_tl, f_yr, f_site:
#' g1 <- make_group_excl_time_state(~ f_tl + f_yr + (1|f_site), df)
#' table(g1)
#'
#' # 2) Drop f_tk and f_tl — group by f_yr and f_site:
#' g2 <- make_group_excl_time_state(~ f_tk + f_tl + f_yr + (1|f_site),
#'                                  df, period = "f_tk", state_var = "f_tl")
#' table(g2)
#'
#' # 3) No remaining variables — returns single-level factor:
#' g3 <- make_group_excl_time_state(~ f_tk + f_tl, df, period = "f_tk", state_var = "f_tl")
#' table(g3)
#'
#' @export
make_group_excl_time_state <- function(fml, data, period = "f_tk", state_var = NULL) {
  bars    <- lme4::findbars(fml)
  re_vars <- if (length(bars)) all.vars(bars[[1]][[3]]) else character(0)

  rhs     <- if (length(fml) == 2) fml[[2]] else fml[[3]]
  fe_vars <- all.vars(rhs)

  exclude_vars <- if (!is.null(state_var)) c(state_var, period) else period
  vars <- setdiff(unique(c(re_vars, fe_vars)), exclude_vars)

  if (length(vars) == 0L) {
    factor(rep(1L, nrow(data)))
  } else {
    interaction(data[vars], drop = TRUE)
  }
}

#' @title Construct Grouping Factor from Fixed and Random Effects
#' @description
#' Combine any factor columns in a model frame (fixed effects) and a list
#' of random-effect factors to produce a single grouping factor via interaction.
#'
#' @param mf_full A `model.frame` containing the fixed-effect factor columns.
#' @param flist A list of factor vectors for random-effect grouping (e.g., from `mkReTrms(...)$flist`).
#'
#' @return A factor of length `nrow(mf_full)` representing the interaction of
#'   all fixed and random grouping factors. If neither fixed nor random factors
#'   are present, returns a single-level factor.
#'
#' @examples
#' # Sample data with two fixed-effect factors
#' df <- data.frame(
#'   f1 = factor(c("A", "A", "B", "B")),
#'   f2 = factor(c("X", "Y", "X", "Y"))
#' )
#' mf_fixed <- model.frame(~ f1 + f2, data = df)
#'
#' # 1) Only fixed effects -> groups by f1 and f2
#' grp1 <- make_process_group(mf_fixed, list())
#' table(grp1)
#'
#' # 2) Only random effects -> ignore fixed, group by random factor
#' rand1 <- factor(c("G1", "G1", "G2", "G2"))
#' grp2 <- make_process_group(mf_fixed[, 0], list(rand1))
#' table(grp2)
#'
#' # 3) Both fixed and random effects -> combine all three
#' grp3 <- make_process_group(mf_fixed, list(rand1))
#' table(grp3)
#'
#' # 4) No factors at all -> single-level grouping
#' empty_mf <- model.frame(~ 1, data = df)
#' grp4 <- make_process_group(empty_mf[, 0], list())
#' table(grp4)
#'
#' @export
make_process_group <- function(mf_full, flist) {
  fixed_fac <- Filter(is.factor, mf_full)
  rand_fac  <- if (length(flist)) as.data.frame(flist) else mf_full[, 0, drop = FALSE]

  if (ncol(fixed_fac) + ncol(rand_fac) == 0) {
    factor(rep(1L, nrow(mf_full)))
  } else {
    interaction(cbind(fixed_fac, rand_fac), drop = TRUE)
  }
}


#' @title Generate Design Matrices and Grouping Structures for Model Processes
#' @description
#' Constructs fixed‐effect design matrices (`X_list`), random‐effect design
#' matrices (`Z_list`), random‐effect term objects (`Z_terms`), grouping
#' indices (`groups`), process‐level grouping factors (`group_factor`),
#' and grouping factors excluding the period variable (`group_noperiod`) for
#' each process in a named list of formulas.
#'
#' @param formulas Named list of formulas, one per process (e.g.,
#'   `list(phi = ~ f_yr + (1|f_site), p = ~ f_yr)`).
#' @param state Not used for exclusion anymore (always drop only the `period` variable).
#' @param data A `data.frame` of covariates, including any variables in `formulas`,
#'   plus the columns required by `make_RTMB_data_list()`.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{group_factor}}{Full joint grouping factors (e.g. sex×state×time).}
#'   \item{\code{group_noperiod}}{Grouping factors dropping only the period (e.g. sex×state).}
#'   \item{…}{(other pieces omitted for brevity).}
#' }
#'
#' @examples
#' #--- common sample data setup:
#' df <- data.frame(
#'   f_sex = factor(rep(c("M","F"), each = 17*3)),
#'   f_tl  = factor(rep(1:3, times = 2*17)),
#'   f_tk  = rep(1:17, times = 2*3),
#'   r_k   = rep(1:17, times = 2*3),
#'   t_l   = f_tl,
#'   f_tl  = f_tl,
#'   r_l   = f_tl,
#'   tag   = seq_len(2*3*17),
#'   n     = sample(1:5, 2*3*17, TRUE)
#' )
#'
#' # -------------------------------------------------------------------
#' # Scenario A: φ random by state AND time  -> 3 sexes×states × 17 periods
#' # -------------------------------------------------------------------
#' formulasA <- list(
#'   phi    = ~ -1 + f_sex + (1|f_tl) + (1|f_tk),
#'   p      = ~ 1,
#'   w      = ~ 1,
#'   t_var  = ~ 1,
#'   Nsuper = ~ 1
#' )
#' designA <- make_design_list(formulasA, state = NULL, data = df)
#' # non-period groups = sex×state  => G_np = 2*3 = 6
#' nlevels(designA$group_noperiod$phi)   # 6
#' # full groups = sex×state×time => G_full = 2*3*17 = 102
#' nlevels(designA$group_factor$phi)     # 102
#'
#' # So your S2 array for Scenario A should be:
#' #   dim = c(G_full=102, n_states, n_states, s=17)
#'
#' # -------------------------------------------------------------------
#' # Scenario B: φ fixed by state, random only by time
#' #    -> keep state in fixed effects, but only random on f_tk
#' # -------------------------------------------------------------------
#' formulasB <- list(
#'   phi    = ~ -1 + f_sex + f_tl + (1|f_tk),
#'   p      = ~ 1,
#'   w      = ~ 1,
#'   t_var  = ~ 1,
#'   Nsuper = ~ 1
#' )
#' designB <- make_design_list(formulasB, state = NULL, data = df)
#' # non-period groups = sex×state  => G_np = 2*3 = 6
#' nlevels(designB$group_noperiod$phi)   # 6
#' # full groups = sex×time (state is fixed so no extra grouping) => 2*17 = 34
#' nlevels(designB$group_factor$phi)     # 34
#'
#' # So your S2 array for Scenario B should be:
#' #   dim = c(G_np=6, n_states, n_states, s=17)
#'
#' @export
make_design_list <- function(formulas, state, data) {
  data <- as.data.frame(data)
  data$._dummy_response <- 0L

  design_list <- list(
    formulas = formulas,
    X_list         = list(),
    Z_list         = list(),
    Z_terms        = list(),
    groups         = list(),
    group_factor   = list(),
    re_factor      = list(),
    group_noperiod = list()
  )

  for (proc in names(formulas)) {
    fml <- formulas[[proc]]
    if (length(fml) == 2)
      fml <- reformulate(deparse(fml[[2]]), response = "._dummy_response")

    mf_full <- model.frame(lme4::subbars(fml), data)
    design_list$X_list[[proc]] <- model.matrix(lme4::nobars(fml), mf_full)

    bars    <- lme4::findbars(fml)
    Zt_trm  <- if (length(bars)) lme4::mkReTrms(bars, mf_full) else NULL

    if (!is.null(Zt_trm)) {
      design_list$Z_terms[[proc]] <- Zt_trm
      design_list$Z_list[[proc]]  <- t(as.matrix(Zt_trm$Zt))
      design_list$groups[[proc]]   <-
        setNames(lapply(Zt_trm$flist, as.integer),
                 paste0("group_", proc, "_re", seq_along(Zt_trm$flist)))
      design_list$re_factor[[proc]] <- factor(rep(1L, nrow(mf_full)))
    } else {
      design_list$Z_terms[[proc]] <- NULL
      design_list$Z_list[[proc]]  <- matrix(0L, nrow(mf_full), 0L)
      design_list$groups[[proc]]   <- list()
      design_list$re_factor[[proc]] <- factor(rep(1L, nrow(mf_full)))
    }

    design_list$group_factor[[proc]] <-
      make_process_group(mf_full, if (is.null(Zt_trm)) list() else Zt_trm$flist)
  }

  for (proc in intersect(names(formulas), c("w", "phi", "p", "t_var", "Nsuper"))) {
    design_list$group_noperiod[[proc]] <-
      make_group_without_period(formulas[[proc]], data, state_var = state)
  }

  design_list
}

#' @title Generate Design Matrices and Grouping Structures (v2)
#' @description
#' Generates X and Z matrices, random effect terms, and grouping factors (including
#' full group, group without period, and group without period and state) for each
#' process in the supplied formulas. Also generates minimal design matrices for each
#' unique group level used to map process parameters into S2/D2 arrays.
#'
#' @param formulas Named list of formulas, one per process.
#' @param input List with names `state` and `input` to identify state and period variables.
#' @param data Data frame containing all covariates used in formulas.
#'
#' @return A list containing components:
#' - formulas: Original formulas
#' - X_obs: Full fixed effect design matrices (for observations)
#' - X_group: Minimal fixed effect design matrices (one per group)
#' - Z_list: Random effect design matrices
#' - Z_terms: Random effect structure
#' - groups: Raw grouping index lists
#' - group_factor: Full grouping factors
#' - group_noperiod: Grouping factors excluding period
#' - group_excl_time_state: Grouping excluding period AND state
#' - data: Input data used
#'
#' @examples
#' df <- data.frame(
#'   f_tk   = rep(1:5, each = 6),
#'   f_tl   = factor(sample(c("A", "B", "C"), 30, replace = TRUE)),
#'   f_yr   = factor(sample(2010:2012, 30, replace = TRUE)),
#'   f_site = factor(sample(1:3, 30, replace = TRUE))
#' )
#'
#' formulas <- list(
#'   phi = ~ f_yr + (1|f_site),
#'   p   = ~ 1
#' )
#'
#' design <- make_design_list_2(formulas, list(state = "f_tl", input = "f_tk"), df)
#' head(design$X_obs$phi)
#' head(design$X_group$phi)
#'
#' @export
make_design_list_1 <- function(formulas, input, data) {
  data <- as.data.frame(data)
  data$._dummy_response <- 0L

  design_list <- list(
    formulas              = formulas,
    X_obs                 = list(),
    X_group               = list(),
    Z_list                = list(),
    Z_terms               = list(),
    groups                = list(),
    group_factor          = list(),
    re_factor             = list(),
    group_noperiod        = list(),
    group_excl_time_state = list(),
    data                  = data
  )

  for (proc in names(formulas)) {
    fml <- formulas[[proc]]
    if (length(fml) == 2)
      fml <- reformulate(deparse(fml[[2]]), response = "._dummy_response")

    mf_full <- model.frame(lme4::subbars(fml), data)
    design_list$X_obs[[proc]] <- model.matrix(lme4::nobars(fml), mf_full)

    bars    <- lme4::findbars(fml)
    Zt_trm  <- if (length(bars)) lme4::mkReTrms(bars, mf_full) else NULL

    if (!is.null(Zt_trm)) {
      design_list$Z_terms[[proc]] <- Zt_trm
      design_list$Z_list[[proc]]  <- t(as.matrix(Zt_trm$Zt))
      design_list$groups[[proc]]  <- setNames(lapply(Zt_trm$flist, as.integer),
                                              paste0("group_", proc, "_re", seq_along(Zt_trm$flist)))
      design_list$re_factor[[proc]] <- factor(rep(1L, nrow(mf_full)))
    } else {
      design_list$Z_terms[[proc]] <- NULL
      design_list$Z_list[[proc]]  <- matrix(0L, nrow(mf_full), 0L)
      design_list$groups[[proc]]  <- list()
      design_list$re_factor[[proc]] <- factor(rep(1L, nrow(mf_full)))
    }

    group_fac <- make_process_group(mf_full, if (is.null(Zt_trm)) list() else Zt_trm$flist)
    design_list$group_factor[[proc]] <- group_fac

    rep_idx <- match(levels(group_fac), group_fac)
    mf_group <- mf_full[rep_idx, , drop = FALSE]
    design_list$X_group[[proc]] <- model.matrix(lme4::nobars(fml), mf_group)
  }

  for (proc in intersect(names(formulas), c("phi", "p", "w", "t_var", "Nsuper"))) {
    design_list$group_noperiod[[proc]] <-
      make_group_without_period(formulas[[proc]], data, state_var = input$state)

    design_list$group_excl_time_state[[proc]] <-
      make_group_without_period(formulas[[proc]], data, period = input$input,
                                state_var = input$state)
  }

  design_list
}

#' @title Map Joint Group Levels to Process-Specific Levels
#' @description
#' For each level of a joint grouping factor, find the corresponding level
#' in a process-specific factor.  Returns an integer vector giving, for each
#' joint level, the index of the first row in which that joint level appears
#' in the process factor.
#'
#' @param joint_fac A factor defining the joint grouping (e.g., interaction of multiple processes).
#' @param proc_fac A factor defining a single process grouping for the same rows.
#'
#' @return An integer vector of length `nlevels(joint_fac)`, where each entry
#'   is the integer code of the corresponding level in `proc_fac`.
#'
#' @examples
#' # Suppose joint grouping has 3 levels, and proc grouping repeats those:
#' joint_fac <- factor(c("A","B","C","A","B","C"))
#' proc_fac  <- factor(c(1,1,2,2,3,3))
#' map_joint_to_proc(joint_fac, proc_fac)
#' #> returns c(1, 1, 2)
#'
#' @export
map_joint_to_proc <- function(joint_fac, proc_fac) {
  sapply(levels(joint_fac), function(lvl) {
    as.integer(proc_fac[ match(TRUE, joint_fac == lvl) ])
  })
}


#' @title Initialize Model Parameters from Design List
#' @description
#' Create starting values for all fixed-effect coefficients (`beta_*`),
#' random effects (`u_*`), and their log-standard deviations (`log_sd_*`)
#' based on the provided design list.  Sensible defaults are used for
#' processes \code{Nsuper}, \code{phi}, and \code{p}, and if \code{data$state}
#' is present, parameters for \code{t_var} are initialized accordingly.
#'
#' @param design A list output by \code{\link{make_design_list}()}, containing
#'   \code{X_list}, \code{Z_list}, \code{group_factor}, etc.
#' @param data The original \code{data.frame} used to build \code{design};
#'   used to detect \code{data$state} for \code{t_var}.
#' @param init_sd Numeric; initial standard deviation for random-effect priors (default: 0.1).
#'
#' @return A named list of starting parameter vectors:
#'   \itemize{
#'     \item \code{beta_<proc>}: zeros for each fixed effect.
#'     \item \code{u_<proc>}: zeros for each random-effect column (or \code{NA}).
#'     \item \code{log_sd_<proc>}: \code{log(init_sd)} where random effects exist (else \code{NA}).
#'   }
#'
#' @examples
#' # Create a dummy design for two processes phi and p:
#' formulas <- list(phi = ~ f_yr + (1|f_site), p = ~ f_yr)
#' df <- data.frame(
#'   f_yr   = factor(rep(2010:2011, each = 4)),
#'   f_site = factor(rep(1:4, times = 2)),
#'   f_tk   = rep(1:4, times = 2),
#'   f_tl   = factor(rep(c("A","B","C","D"), times = 2)),
#'   r_k    = rep(1:4, times = 2),
#'   r_l    = factor(rep(1:4, times = 2)),
#'   tag    = 1:8,
#'   n      = sample(1:5, 8, replace = TRUE)
#' )
#' design <- make_design_list(formulas, state = "f_tl", data = df)
#' params <- make_all_params(design, df, init_sd = 0.2)
#' names(params)
#'
#' @export
make_all_params <- function(design, data, init_sd = 0.1) {
  params <- list()

  for (proc in names(design$group_factor)) {
    X <- design$X_list[[proc]]
    Z <- design$Z_list[[proc]]

    params[[paste0("beta_", proc)]] <-
      if (ncol(X) > 0) rep(0, ncol(X)) else NA

    if (ncol(Z) > 0) {
      params[[paste0("u_", proc)]]      <- rep(0, ncol(Z))
      params[[paste0("log_sd_", proc)]] <- log(init_sd)
    } else {
      params[[paste0("u_", proc)]]      <- NA
      params[[paste0("log_sd_", proc)]] <- NA
    }

    # For t_var when data$state exists, override defaults
    if (!is.null(data$state) && proc == "t_var") {
      params[[paste0("u_", proc)]]      <- 0
      params[[paste0("log_sd_", proc)]] <- 1
    }
  }

  # Sensible defaults for intercepts
  if ("beta_Nsuper" %in% names(params)) {
    params$beta_Nsuper <- rep(10, length(params$beta_Nsuper))
  }
  if ("beta_phi" %in% names(params)) {
    params$beta_phi    <- rep(-2, length(params$beta_phi))
  }
  if ("beta_p" %in% names(params)) {
    params$beta_p      <- rep(-1, length(params$beta_p))
  }

  # If a state dimension is used, set up t_var random effects
  if (!is.null(data$state)) {
    n_states <- 5
    n_par_M  <- n_states * (n_states - 1) / 2
    params$u_t_var <- rep(0, n_par_M)
  }

  params
}


#' @title Construct RTMB Data List from Design and Capture Data
#' @description
#' Centralizes construction of the full data list required by an RTMB model.
#' Pulls together observed capture histories, design matrices, grouping factors,
#' joint-to-process mappings, and auxiliary summaries (e.g., counts per time,
#' spatial proportions) into a single named list.
#'
#' @param design A design list generated by \code{\link{make_design_list}()},
#'   containing \code{X_list}, \code{Z_list}, \code{group_factor}, and
#'   \code{group_noperiod} components.
#' @param data A \code{data.frame} of individual‐level capture data.  Must
#'   include columns \code{t_k}, \code{r_k}, \code{t_l}, \code{f_tl},
#'   \code{r_l}, \code{tag}, and \code{n}, plus any covariates used in
#'   \code{design}.
#' @param state Optional string naming the state column (e.g., \code{"t_l"}).
#'   If provided, \code{state} is stored in the output and used for spatial
#'   summaries.
#' @param period Character. Column name for the period/time variable
#'   (default: \code{"t_k"}); used to tabulate counts per time step.
#'
#' @return A named list with components:
#' \describe{
#'   \item{\code{t_k}, \code{r_k}, \code{t_l}, \code{f_tl}, \code{r_l}, \code{tag}, \code{n}}{Raw data vectors.}
#'   \item{\code{s}}{Maximum observed time or recapture index.}
#'   \item{\code{period_fac}}{Factor of the period variable.}
#'   \item{\code{X_list}, \code{Z_list}}{Design matrices from \code{design}.}
#'   \item{\code{group_factor}}{List of joint grouping factors.}
#'   \item{\code{phi_for_joint}, \code{p_for_joint}, \code{w_for_joint}, \code{tvar_for_joint}, \code{Nsuper_for_joint}}{Integer mappings from joint to process groups.}
#'   \item{\code{joint_group}, \code{Nsuper_group}}{Factors of joint and Nsuper groups.}
#'   \item{\code{G}, \code{K}}{Number of joint groups and number of Nsuper levels.}
#'   \item{\code{A}}{Membership matrix (K × G) mapping Nsuper levels to joint groups.}
#'   \item{\code{u}}{Vector of counts tagged per period.}
#'   \item{\code{uTot_g}}{Total counts per Nsuper group.}
#'   \item{\code{sp_pr}}{Spatial proportions of detections.}
#'   \item{\code{state}}{Same as input \code{state}.}
#' }
#'
#' @examples
#' # 1) Build a simple design for processes phi and p
#' formulas <- list(
#'   phi = ~ f_yr + (1|f_site),
#'   p   = ~ f_yr
#' )
#' df <- data.frame(
#'   f_yr   = factor(rep(2010:2011, each = 4)),
#'   f_site = factor(rep(1:4, times = 2)),
#'   t_k    = rep(1:4, times = 2),
#'   r_k    = rep(1:4, times = 2),
#'   t_l    = factor(rep(c("A","B","C","D"), times = 2)),
#'   f_tl   = factor(rep(c("A","B","C","D"), times = 2)),
#'   r_l    = factor(rep(1:4, times = 2)),
#'   tag    = 1:8,
#'   n      = sample(1:5, 8, replace = TRUE)
#' )
#'
#' # 2) Create design list, dropping f_tl for group_noperiod
#' design <- make_design_list(formulas, state = "t_l", data = df)
#'
#' # 3) Build RTMB data list
#' rtmb_data <- make_RTMB_data_list(design, df, state = "t_l", period = "t_k")
#'
#' # Inspect key elements:
#' names(rtmb_data)
#' rtmb_data$G
#' head(rtmb_data$u)
#' rtmb_data$sp_pr
#'
#' @export
make_RTMB_data_list <- function(design, data, state = NULL, period = "t_k") {

  phi_group     <- design$group_noperiod[['phi']]
  p_group       <- design$group_noperiod[['p']]
  w_group       <- design$group_noperiod[['w']]
  tvar_group    <- design$group_noperiod[['t_var']]
  Nsuper_group  <- design$group_noperiod[['Nsuper']]

  joint_group <- interaction(phi_group, p_group, w_group, tvar_group, drop = TRUE)
  G           <- nlevels(joint_group)

  phi_for_joint      <- map_joint_to_proc(joint_group, phi_group)
  p_for_joint        <- map_joint_to_proc(joint_group, p_group)
  w_for_joint        <- map_joint_to_proc(joint_group, w_group)
  tvar_for_joint     <- map_joint_to_proc(joint_group, tvar_group)
  Nsuper_for_joint   <- map_joint_to_proc(joint_group, Nsuper_group)

  if (nlevels(Nsuper_group) == 1) {
    K <- 1L
    A <- matrix(1, nrow = 1, ncol = G)
  } else {
    K <- nlevels(Nsuper_group)
    A <- Matrix::sparseMatrix(
      i    = as.integer(Nsuper_group),
      j    = as.integer(joint_group),
      x    = 1,
      dims = c(K, G)
    )
    A <- as.matrix(A)
  }

  # Counts per period
  s     <- max(na.omit(c(data[[period]], data$r_k)))
  u_vec <- numeric(s)
  u_sum <- aggregate(list(n = data$n), by = list(t = data[[period]]), sum)
  u_vec[u_sum$t] <- u_sum$n

  # Spatial proportions
  sp_count <- aggregate(list(n = data$n),
                        by = list(loc = data[[ state %||% 't_l' ]]), sum)
  sp_pr <- sp_count$n / sum(sp_count$n)

  # Totals per Nsuper group
  if (length(Nsuper_group) > 1) {
    uTot_g <- as.numeric(tapply(data$n, Nsuper_group, sum))
  } else {
    uTot_g <- sum(data$n)
  }

  list(
    design = design,
    t_k               = data$t_k,
    r_k               = data$r_k,
    t_l               = data$t_l,
    f_tl              = data$f_tl,
    r_l               = data$r_l,
    tag               = as.integer(data$tag),
    n                 = data$n,
    s                 = s,
    period_fac        = data[[period]],
    X_list            = design$X_list,
    Z_list            = design$Z_list,
    group_factor      = design$group_factor,
    phi_for_joint     = phi_for_joint,
    p_for_joint       = p_for_joint,
    w_for_joint       = w_for_joint,
    tvar_for_joint    = tvar_for_joint,
    Nsuper_for_joint  = Nsuper_for_joint,
    joint_group       = joint_group,
    Nsuper_group      = Nsuper_group,
    G                 = G,
    K                 = K,
    A                 = A,
    u                 = u_vec,
    uTot_g            = uTot_g,
    sp_pr             = sp_pr,
    state             = state
  )
}

#' @title Create Data Frame of Joint-to-Process Group Mappings
#' @description
#' Generates a tidy \code{data.frame} that shows, for each level of a joint
#' grouping factor, the corresponding level in one or more process-specific
#' grouping factors.
#'
#' @param joint_group A factor representing the joint grouping (e.g., interaction of processes).
#' @param group_factor A named list of factor vectors (one per process) aligning
#'   with the same rows as \code{joint_group}.
#'
#' @return A \code{data.frame} with one row per joint level and one column
#'   named \code{joint}, plus one column per element of \code{group_factor},
#'   giving the integer code of the process-specific group for that joint level.
#'
#' @examples
#' joint_group <- factor(c("A","B","A","C","B"))
#' group_factor <- list(
#'   phi = factor(c(1,2,1,3,2)),
#'   p   = factor(c(1,1,1,2,2))
#' )
#' make_group_mapping_df(joint_group, group_factor)
#'
#' @export
make_group_mapping_df <- function(joint_group, group_factor) {
  G <- nlevels(joint_group)
  data.frame(
    joint = levels(joint_group),
    lapply(group_factor, function(gf) {
      vapply(levels(joint_group), function(lvl) {
        idx <- which(joint_group == lvl)[1]
        as.integer(gf[idx])
      }, integer(1))
    })
  )
}


#' @title Plot Sankey-Style Mapping of Joint to Process Groups
#' @description
#' Visualizes how joint grouping levels map to each process-specific grouping
#' by drawing a tile grid colored by group and labeled with group indices.
#'
#' @param data A list containing at least these elements:
#'   \describe{
#'     \item{\code{joint_group}}{factor of joint groups}
#'     \item{\code{phi_for_joint}, \code{p_for_joint}, \code{w_for_joint},
#'           \code{tvar_for_joint}, \code{Nsuper_for_joint}}{integer vectors mapping each joint level to a process group}
#'   }
#'
#' @return A \code{ggplot} object displaying the group mapping.
#'
#' @examples
#' # Suppose rtmb_data is returned by make_RTMB_data_list()
#' # and contains the necessary mapping vectors:
#' # plot_Sankey(rtmb_data)
#'
#' @export
plot_Sankey <- function(data) {
  library(viridis)
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(scales)

  mapping_df <- data.frame(
    joint   = levels(data$joint_group),
    phi     = data$phi_for_joint,
    p       = data$p_for_joint,
    w       = data$w_for_joint,
    t_var   = data$tvar_for_joint,
    Nsuper  = data$Nsuper_for_joint
  ) %>%
    pivot_longer(cols = -joint, names_to = "process", values_to = "group") %>%
    mutate(
      joint   = factor(joint, levels = unique(joint)),
      process = factor(process, levels = c("phi", "p", "w", "t_var", "Nsuper")),
      group   = as.integer(group)
    ) %>%
    group_by(process) %>%
    mutate(group_scaled = rescale(group, to = c(1, 100))) %>%
    ungroup()

  ggplot(mapping_df, aes(x = process, y = joint)) +
    geom_tile(aes(fill = group_scaled), color = "white") +
    geom_text(aes(label = group), size = 2.5) +
    scale_fill_viridis_c(name = "Group", option = "plasma") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 7),
      panel.grid = element_blank()
    ) +
    labs(
      title = "Group Mapping by Process",
      x     = "Process",
      y     = "Joint Group"
    )
}


#' @title Identify Random-Effect Parameter Names
#' @description
#' Examines a named parameter list and returns the names of those elements
#' corresponding to random-effect vectors (length > 1).
#'
#' @param params A named list of model parameters (as from \code{\link{make_all_params}()}).
#'
#' @return A character vector of parameter names that are random-effect vectors.
#'
#' @examples
#' params <- list(
#'   u_phi = rnorm(3),
#'   u_p   = numeric(0),
#'   beta_phi = 0,
#'   u_w   = rnorm(2)
#' )
#' make_random(params)
#'
#' @export
make_random <- function(params) {
  re_names <- c("u_phi", "u_p", "u_w", "u_v", "u_Nsuper", "u_t_var")
  names(Filter(function(x) length(x) > 1, params[re_names]))
}


#' @title Construct Parameter Map for TMB
#' @description
#' Builds a \code{map} list indicating which parameters should be held fixed
#' (mapped to \code{NA}) versus estimated, based on the provided parameter list.
#'
#' @param params A named list of parameters (from \code{\link{make_all_params}()}).
#'
#' @return A named list of \code{factor} objects suitable for passing as the
#'   \code{map} argument in \code{TMB::MakeADFun()}.
#'
#' @examples
#' design <- make_design_list(list(phi = ~1), state = NULL, data = data.frame(t_k=1,r_k=1,t_l=1,f_tl=1,r_l=1,tag=1,n=1))
#' params <- make_all_params(design, data.frame(state=NULL))
#' make_map(params)
#'
#' @export
make_map <- function(params) {
  map <- list()

  # Always fixed entries (customize as needed)
  map$beta_v       <- factor(NA)
  map$beta_w       <- factor(NA)
  map$u_v          <- factor(NA)
  map$log_sd_v     <- factor(NA)
  map$u_Nsuper     <- factor(NA)
  map$log_sd_Nsuper<- factor(NA)
  map$beta_t_var   <- factor(NA)

  if (is.na(params$log_sd_t_var)) {
    map$u_t_var     <- factor(rep(NA, length(params$u_t_var)))
    map$log_sd_t_var<- factor(NA)
  }
  if (is.na(params$log_sd_phi)) {
    map$u_phi       <- factor(rep(NA, length(params$u_phi)))
    map$log_sd_phi  <- factor(NA)
  }
  if (is.na(params$log_sd_p)) {
    map$u_p         <- factor(rep(NA, length(params$u_p)))
    map$log_sd_p    <- factor(NA)
  }
  if (is.na(params$log_sd_w)) {
    map$u_w         <- factor(rep(NA, length(params$u_w)))
    map$log_sd_w    <- factor(NA)
  }

  map
}

#’ @title Map Each Group to Its Observed Locations
#’ @description
#’ Given a grouping factor (length = nrows of your data) and a vector of
#’ location codes (same length), returns a named list where each entry is
#’ the unique set of locations in which that group was observed.
#’
#’ @param group_fac A factor of length N (e.g. `design$group_noperiod$phi`).
#’ @param loc_var   A vector (numeric or factor) of length N with location IDs.
#’
#’ @return A named list of length `nlevels(group_fac)`; each element is a
#’   vector of the unique `loc_var` values for that group.
#’
#’ @examples
#’ # sample data
#’ df <- data.frame(
#’   sex   = factor(c("F","F","M","M","J","J","J")),
#’   state = factor(c(1,2,1,3,2,3,4))
#’ )
#’ # pretend group_noperiod = sex
#’ grp <- df$sex
#’ # locations of each tag in df$state
#’ lookup <- map_group_to_locs(grp, df$state)
#’ lookup
#’ #> $F
#’ #> [1] 1 2
#’ #> $M
#’ #> [1] 1 3
#’ #> $J
#’ #> [1] 2 3 4
#’
#’ @export
map_group_to_locs <- function(group_fac, loc_var) {
  stats::setNames(
    lapply(levels(group_fac), function(lvl) {
      unique(loc_var[ group_fac == lvl ])
    }),
    levels(group_fac)
  )
}

#' @title Build a Full Design Grid for Any φ‐Formula
#' @description
#' Given a φ‐formula that may include any number of fixed & random effects,
#' plus a state variable (`f_tl`) and a period variable (`f_tk`),
#' build a data.frame with one row per possible combination of:
#'   - every factor in the formula (fixed or grouping),
#'   - the full set of state levels, and
#'   - the full set of period values.
#'
#' @param phi_formula A one‐sided formula for φ, e.g. `~1 + f_sex + (1|f_tl)`.
#' @param data        The original data frame (must contain all relevant columns).
#' @param state       Optional name of the state variable (e.g. `"f_tl"`).
#' @param period      Name of the period variable (e.g. `"f_tk"`).
#'
#' @return A data.frame whose columns include every factor in `phi_formula`,
#'   plus `state` and `period`, with one row per combination of their levels.
#'
#' @examples
#' # Suppose your data has factors f_sex (3 levels), f_tl (5 levels), f_tk in 1:17
#' full_grid <- create_full_grid(~ 1 + f_sex + (1|f_tl), d,
#'                              state = "f_tl", period = "f_tk")
#' dim(full_grid)
#' # 3 × 5 × 17 = 255 rows
#'
#' # You can now do:
#' fe_phi <- lme4::nobars(~ 1 + f_sex + (1|f_tl))
#' X_full  <- model.matrix(fe_phi, full_grid)
#' bars    <- lme4::findbars(~ 1 + f_sex + (1|f_tl))
#' Z_full  <- t(as.matrix(lme4::mkReTrms(bars, full_grid)$Zt))
#' # X_full is 255×p,   Z_full is 255×q
#'
#' @export
create_full_grid <- function(phi_formula, data, state = NULL, period = "f_tk") {
  # 1) fixed‐effect vars
  fe <- lme4::nobars(phi_formula)
  fe_vars <- all.vars(fe[[length(fe)]])
  # 2) random‐effect grouping vars
  bars   <- lme4::findbars(phi_formula)
  re_vars <- unique(unlist(lapply(bars, function(b) all.vars(b[[3]]))))
  # 3) union & drop period/state if present
  vars   <- setdiff(unique(c(fe_vars, re_vars)), c(period, state))

  # 4) keep only those that are factors in your data
  facs   <- vars[sapply(vars, function(v) v %in% names(data) && is.factor(data[[v]]))]
  grid_list <- setNames(
    lapply(facs, function(v) levels(data[[v]])),
    facs
  )

  # 5) add state & period back in
  if (!is.null(state)) {
    if (!state %in% names(data) || !is.factor(data[[state]]))
      stop("`state` must be a factor column in data")
    grid_list[[state]] <- levels(data[[state]])
  }
  if (!period %in% facs) {
    if (!period %in% names(data))
      stop("`period` not found in data")
    # treat period as integer 1..max
    if (is.factor(data[[period]])) {
      grid_list[[period]] <- levels(data[[period]])
    } else {
      grid_list[[period]] <- seq_len(max(data[[period]], na.rm = TRUE))
    }
  }

  # 6) expand
  expand.grid(grid_list, stringsAsFactors = TRUE)
}

#' @title Sanity Check for Process Lookup and Parameter Mapping
#' @description
#' Validates that the lookup table for a process (e.g., `phi` or `p`) is consistent
#' with the expected dimensions of the S2 or D2 arrays and that all index combinations
#' are assigned appropriate parameter rows from the minimal design matrix.
#'
#' This function is useful for validating your design matrix and parameter map setup
#' before embedding the logic in RTMB.
#'
#' @param proc Character. Process name (e.g., "phi").
#' @param design The design list from `make_design_list_2()`.
#' @param lookup The output from `make_process_lookup()`.
#' @param params A list of parameter vectors from `make_param_vectors()`.
#' @param state_var Character. Name of the state variable.
#' @param period_var Character. Name of the period variable (default = "f_tk").
#'
#' @return A list containing:
#' - `dim_S`: Dimensions of the process matrix (G x S x S x T)
#' - `design_rows`: Table of how many times each row in the minimal design matrix is used
#' - `valid`: Logical flag indicating if all rows referenced are within the bounds of `params`
#' - `table`: A summary table of usage counts by design row
#'
#' @examples
#' set.seed(123)
#' df <- data.frame(
#'   f_tk   = rep(1:5, each = 12),
#'   f_tl   = factor(rep(c("A", "B", "C"), times = 20)),
#'   f_yr   = factor(sample(rep(2010:2012, length.out = 60))),
#'   f_site = factor(sample(rep(1:3, length.out = 60)))
#' )
#'
#' formulas <- list(
#'   phi = ~ f_yr + (1|f_site),
#'   p   = ~ 1
#' )
#'
#' design <- make_design_list_2(formulas, list(state = "f_tl", input = "f_tk"), data = df)
#' lookup_phi <- make_process_lookup("phi", design, df, state_var = "f_tl",
#'                                   full_groups = design$group_factor$phi)
#' params_phi <- make_param_vectors("phi", design)
#' sanity_check_process("phi", design, lookup_phi, params_phi, state_var = "f_tl")
#'
#' @export
sanity_check_process <- function(proc, design, lookup, params,
                                 state_var = NULL, period_var = "f_tk") {
  # Determine sizes
  G <- max(lookup$g)
  S <- if (!is.null(state_var)) length(levels(design$data[[state_var]])) + 1L else 2L
  T <- max(lookup$time)

  design_rows <- lookup$design_row
  beta_len <- length(params$beta)

  # Check if all referenced rows are within bounds
  valid_refs <- all(design_rows >= 1 & design_rows <= beta_len)

  usage_table <- table(design_rows)

  list(
    dim_S = c(G, S, S, T),
    design_rows = design_rows,
    valid = valid_refs,
    table = usage_table
  )
}

#' @title Generate Design Matrices and Grouping Structures (v2)
#' @description
#' Generates X and Z matrices, random effect terms, and grouping factors (including
#' full group, group without period, and group without period and state) for each
#' process in the supplied formulas. Also generates minimal design matrices for each
#' unique group level used to map process parameters into S2/D2 arrays.
#'
#' @param formulas Named list of formulas, one per process.
#' @param input List with names `state` and `input` to identify state and period variables.
#' @param data Data frame containing all covariates used in formulas.
#'
#' @return A list containing components:
#' - formulas: Original formulas
#' - X_obs: Full fixed effect design matrices (for observations)
#' - X_group: Minimal fixed effect design matrices (one per group)
#' - Z_obs: Full random effect design matrices
#' - Z_group: Minimal random effect design matrices (one per group)
#' - group_factor: Full grouping factors
#' - group_noperiod: Grouping factors excluding period
#' - group_excl_time_state: Grouping excluding period AND state
#' - data: The input data used
#'
#' @examples
#' df <- data.frame(
#'   f_tk   = rep(1:5, each = 6),
#'   f_tl   = factor(sample(c("A", "B", "C"), 30, replace = TRUE)),
#'   f_yr   = factor(sample(2010:2012, 30, replace = TRUE)),
#'   f_site = factor(sample(1:3, 30, replace = TRUE))
#' )
#'
#' formulas <- list(
#'   phi = ~ f_yr + (1|f_site),
#'   p   = ~ 1
#' )
#'
#' design <- make_design_list_2(formulas, list(state = "f_tl", input = "f_tk"), df)
#' head(design$X_obs$phi)
#' head(design$X_group$phi)
#' head(design$Z_obs$phi)
#' head(design$Z_group$phi)
#'
#' @export
make_design_list_2 <- function(formulas, input, data) {
  data <- as.data.frame(data)
  data$._dummy_response <- 0L

  design_list <- list(
    formulas              = formulas,
    X_obs                 = list(),
    X_group               = list(),
    Z_obs                 = list(),
    Z_group               = list(),
    Z_terms               = list(),
    groups                = list(),
    group_factor          = list(),
    re_factor             = list(),
    group_noperiod        = list(),
    group_excl_time_state = list(),
    data                  = data,
    state                 = input$state
  )

  for (proc in names(formulas)) {
    fml <- formulas[[proc]]
    if (length(fml) == 2)
      fml <- reformulate(deparse(fml[[2]]), response = "._dummy_response")

    mf_full <- model.frame(lme4::subbars(fml), data)
    design_list$X_obs[[proc]] <- model.matrix(lme4::nobars(fml), mf_full)

    bars    <- lme4::findbars(fml)
    Zt_trm  <- if (length(bars)) lme4::mkReTrms(bars, mf_full) else NULL

    if (!is.null(Zt_trm)) {
      design_list$Z_terms[[proc]] <- Zt_trm
      design_list$Z_obs[[proc]]   <- t(as.matrix(Zt_trm$Zt))
      design_list$groups[[proc]]  <- setNames(lapply(Zt_trm$flist, as.integer),
                                              paste0("group_", proc, "_re", seq_along(Zt_trm$flist)))
      design_list$re_factor[[proc]] <- factor(rep(1L, nrow(mf_full)))
    } else {
      design_list$Z_terms[[proc]] <- NULL
      design_list$Z_obs[[proc]]   <- matrix(0L, nrow(mf_full), 0L)
      design_list$groups[[proc]]  <- list()
      design_list$re_factor[[proc]] <- factor(rep(1L, nrow(mf_full)))
    }

    group_fac <- make_process_group(mf_full, if (is.null(Zt_trm)) list() else Zt_trm$flist)
    design_list$group_factor[[proc]] <- group_fac

    rep_idx <- match(levels(group_fac), group_fac)
    mf_group <- mf_full[rep_idx, , drop = FALSE]
    design_list$X_group[[proc]] <- model.matrix(lme4::nobars(fml), mf_group)

    if (!is.null(Zt_trm)) {
      Zt_trm_group <- lme4::mkReTrms(bars, mf_group)
      design_list$Z_group[[proc]] <- t(as.matrix(Zt_trm_group$Zt))
    } else {
      design_list$Z_group[[proc]] <- matrix(0L, nrow(mf_group), 0L)
    }
  }

  for (proc in intersect(names(formulas), c("phi", "p", "w", "t_var", "Nsuper"))) {
    design_list$group_noperiod[[proc]] <-
      make_group_without_period(formulas[[proc]], data, state_var = input$state)

    design_list$group_excl_time_state[[proc]] <-
      make_group_without_period(formulas[[proc]], data, period = input$input,
                                state_var = input$state)
  }

  design_list
}
