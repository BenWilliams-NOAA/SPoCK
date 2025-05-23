% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Get_Reference_Points.R
\name{Get_Reference_Points}
\alias{Get_Reference_Points}
\title{Wrapper function to get reference points}
\usage{
Get_Reference_Points(
  data,
  rep,
  SPR_x = NULL,
  t_spwn = 0,
  type = "single_region",
  what = "SPR"
)
}
\arguments{
\item{data}{Data list from RTMB}

\item{rep}{Report list from RTMB}

\item{SPR_x}{SPR percentage to target}

\item{t_spwn}{specified mortality time up until spawning}

\item{type}{Whether this is a "single_region" reference point (options are being developed)}

\item{what}{What kind of reference point "SPR" (options are being developed)}
}
\value{
A RTMB list object with report information on estiamted reference points
}
\description{
Wrapper function to get reference points
}
\examples{
\dontrun{
f_40 <- Get_Reference_Points(data = data,
rep = rep,
SPR_x = 0.4,
t_spwn = 0,
type = "single_region",
what = "SPR")
f_40$rep$F_x # estimated reference point
}
}
