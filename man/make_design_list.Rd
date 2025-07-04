% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utlis.R
\name{make_design_list}
\alias{make_design_list}
\title{Generate Design Matrices and Grouping Structures for Model Processes}
\usage{
make_design_list(formulas, state, data)
}
\arguments{
\item{formulas}{Named list of formulas, one per process (e.g.,
`list(phi = ~ f_yr + (1|f_site), p = ~ f_yr)`).}

\item{state}{Not used for exclusion anymore (always drop only the `period` variable).}

\item{data}{A `data.frame` of covariates, including any variables in `formulas`,
plus the columns required by `make_RTMB_data_list()`.}
}
\value{
A list with components:
\describe{
  \item{\code{group_factor}}{Full joint grouping factors (e.g. sex×state×time).}
  \item{\code{group_noperiod}}{Grouping factors dropping only the period (e.g. sex×state).}
  \item{…}{(other pieces omitted for brevity).}
}
}
\description{
Constructs fixed‐effect design matrices (`X_list`), random‐effect design
matrices (`Z_list`), random‐effect term objects (`Z_terms`), grouping
indices (`groups`), process‐level grouping factors (`group_factor`),
and grouping factors excluding the period variable (`group_noperiod`) for
each process in a named list of formulas.
}
\examples{
#--- common sample data setup:
df <- data.frame(
  f_sex = factor(rep(c("M","F"), each = 17*3)),
  f_tl  = factor(rep(1:3, times = 2*17)),
  f_tk  = rep(1:17, times = 2*3),
  r_k   = rep(1:17, times = 2*3),
  t_l   = f_tl,
  f_tl  = f_tl,
  r_l   = f_tl,
  tag   = seq_len(2*3*17),
  n     = sample(1:5, 2*3*17, TRUE)
)

# -------------------------------------------------------------------
# Scenario A: φ random by state AND time  -> 3 sexes×states × 17 periods
# -------------------------------------------------------------------
formulasA <- list(
  phi    = ~ -1 + f_sex + (1|f_tl) + (1|f_tk),
  p      = ~ 1,
  w      = ~ 1,
  t_var  = ~ 1,
  Nsuper = ~ 1
)
designA <- make_design_list(formulasA, state = NULL, data = df)
# non-period groups = sex×state  => G_np = 2*3 = 6
nlevels(designA$group_noperiod$phi)   # 6
# full groups = sex×state×time => G_full = 2*3*17 = 102
nlevels(designA$group_factor$phi)     # 102

# So your S2 array for Scenario A should be:
#   dim = c(G_full=102, n_states, n_states, s=17)

# -------------------------------------------------------------------
# Scenario B: φ fixed by state, random only by time
#    -> keep state in fixed effects, but only random on f_tk
# -------------------------------------------------------------------
formulasB <- list(
  phi    = ~ -1 + f_sex + f_tl + (1|f_tk),
  p      = ~ 1,
  w      = ~ 1,
  t_var  = ~ 1,
  Nsuper = ~ 1
)
designB <- make_design_list(formulasB, state = NULL, data = df)
# non-period groups = sex×state  => G_np = 2*3 = 6
nlevels(designB$group_noperiod$phi)   # 6
# full groups = sex×time (state is fixed so no extra grouping) => 2*17 = 34
nlevels(designB$group_factor$phi)     # 34

# So your S2 array for Scenario B should be:
#   dim = c(G_np=6, n_states, n_states, s=17)

}
