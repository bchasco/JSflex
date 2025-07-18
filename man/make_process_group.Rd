% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utlis.R
\name{make_process_group}
\alias{make_process_group}
\title{Construct Grouping Factor from Fixed and Random Effects}
\usage{
make_process_group(mf_full, flist)
}
\arguments{
\item{mf_full}{A `model.frame` containing the fixed-effect factor columns.}

\item{flist}{A list of factor vectors for random-effect grouping (e.g., from `mkReTrms(...)$flist`).}
}
\value{
A factor of length `nrow(mf_full)` representing the interaction of
  all fixed and random grouping factors. If neither fixed nor random factors
  are present, returns a single-level factor.
}
\description{
Combine any factor columns in a model frame (fixed effects) and a list
of random-effect factors to produce a single grouping factor via interaction.
}
\examples{
# Sample data with two fixed-effect factors
df <- data.frame(
  f1 = factor(c("A", "A", "B", "B")),
  f2 = factor(c("X", "Y", "X", "Y"))
)
mf_fixed <- model.frame(~ f1 + f2, data = df)

# 1) Only fixed effects -> groups by f1 and f2
grp1 <- make_process_group(mf_fixed, list())
table(grp1)

# 2) Only random effects -> ignore fixed, group by random factor
rand1 <- factor(c("G1", "G1", "G2", "G2"))
grp2 <- make_process_group(mf_fixed[, 0], list(rand1))
table(grp2)

# 3) Both fixed and random effects -> combine all three
grp3 <- make_process_group(mf_fixed, list(rand1))
table(grp3)

# 4) No factors at all -> single-level grouping
empty_mf <- model.frame(~ 1, data = df)
grp4 <- make_process_group(empty_mf[, 0], list())
table(grp4)

}
