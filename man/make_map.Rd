% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utlis.R
\name{make_map}
\alias{make_map}
\title{Construct Parameter Map for TMB}
\usage{
make_map(params)
}
\arguments{
\item{params}{A named list of parameters (from \code{\link{make_all_params}()}).}
}
\value{
A named list of \code{factor} objects suitable for passing as the
  \code{map} argument in \code{TMB::MakeADFun()}.
}
\description{
Builds a \code{map} list indicating which parameters should be held fixed
(mapped to \code{NA}) versus estimated, based on the provided parameter list.
}
\examples{
design <- make_design_list(list(phi = ~1), state = NULL, data = data.frame(t_k=1,r_k=1,t_l=1,f_tl=1,r_l=1,tag=1,n=1))
params <- make_all_params(design, data.frame(state=NULL))
make_map(params)

}
