% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utlis.R
\name{prep_param_input}
\alias{prep_param_input}
\title{Prepare Group-Based Predictors and Group-to-State Maps for Model Processes}
\usage{
prep_param_input(
  proc,
  data,
  design,
  beta_list,
  u_list,
  state_var = NULL,
  period_var = "f_tk"
)
}
\arguments{
\item{proc}{String name of the model process (`"phi"`, `"p"`, etc.)}

\item{data}{A `data.frame` of covariates, with time and state columns as needed.}

\item{design}{Output of `make_design_list_2()` (includes formulas and group factors).}

\item{beta_list}{Named list of fixed effect parameter vectors (`beta$phi`, etc.)}

\item{u_list}{Named list of random effect parameter vectors (`u$phi`, etc.)}

\item{state_var}{Optional string naming the state variable (e.g., `"f_tl"`); if `NULL`,
all individuals are assumed to belong to state 1.}

\item{period_var}{String naming the period variable (default `"f_tk"`).}
}
\value{
A list with the following components:
\describe{
  \item{\code{eta}}{Numeric vector of linear predictors for each group.}
  \item{\code{g_to_locs}}{List mapping each group to one or more location indices.}
  \item{\code{group_fac}}{Factor used to define the unique grouping structure.}
  \item{\code{n_groups}}{Number of unique groups.}
  \item{\code{mapping}}{A data.frame of all group × state × period combinations.}
}
}
\description{
General utility for extracting linear predictors and group-to-location mappings
for any model process (`phi`, `p`, `w`, etc.). It computes:
- the linear predictor (`eta`) for each unique group
- the mapping from each group to one or more states (locations)
- a complete index grid over group × state × period combinations

The function ensures that:
- the **period variable** (e.g., `f_tk`) never contributes to the grouping structure
- the **state variable** (e.g., `f_tl`) is excluded from grouping if provided,
  though still used for group-to-location mapping

This precomputes the structure required to fully populate matrices like `S2` and `D2`
even in the presence of missing group × state × period combinations.
}
\examples{
# Simulated dataset
set.seed(42)
df <- data.frame(
  f_tk   = rep(1:5, times = 6),
  f_tl   = factor(rep(1:2, each = 15)),
  f_yr   = factor(rep(2010:2012, length.out = 30)),
  f_site = factor(sample(1:3, 30, TRUE))
)

# Example formulas
formulas <- list(
  phi = ~ f_yr + (1|f_site),
  p   = ~ 1
)

design <- make_design_list_2(formulas, list(state = "f_tl", input = "f_tk"), data = df)
beta  <- list(phi = rep(0.1, ncol(design$X_list$phi)), p = 0)
u     <- list(phi = rep(0, ncol(design$Z_list$phi)),   p = numeric(0))

param_phi <- prep_param_input("phi", df, design, beta, u, state_var = "f_tl")
param_phi$eta         # linear predictor for each group
param_phi$g_to_locs   # locations each group covers
param_phi$mapping     # all group × state × time combinations

}
