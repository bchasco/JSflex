% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utlis.R
\name{plot_Sankey}
\alias{plot_Sankey}
\title{Plot Sankey-Style Mapping of Joint to Process Groups}
\usage{
plot_Sankey(data)
}
\arguments{
\item{data}{A list containing at least these elements:
\describe{
  \item{\code{joint_group}}{factor of joint groups}
  \item{\code{phi_for_joint}, \code{p_for_joint}, \code{w_for_joint},
        \code{tvar_for_joint}, \code{Nsuper_for_joint}}{integer vectors mapping each joint level to a process group}
}}
}
\value{
A \code{ggplot} object displaying the group mapping.
}
\description{
Visualizes how joint grouping levels map to each process-specific grouping
by drawing a tile grid colored by group and labeled with group indices.
}
\examples{
# Suppose rtmb_data is returned by make_RTMB_data_list()
# and contains the necessary mapping vectors:
# plot_Sankey(rtmb_data)

}
