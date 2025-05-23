% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Setup_Fishery.R
\name{Setup_Mod_Catch_and_F}
\alias{Setup_Mod_Catch_and_F}
\title{Setup fishing mortality and catch observations}
\usage{
Setup_Mod_Catch_and_F(
  input_list,
  ObsCatch,
  Catch_Type,
  UseCatch,
  Use_F_pen = 1,
  est_all_regional_F = 1,
  Catch_Constant = NULL,
  sigmaC_spec = "fix",
  sigmaF_spec = "fix",
  sigmaF_agg_spec = "fix",
  ...
)
}
\arguments{
\item{input_list}{List containing a data list, parameter list, and map list}

\item{ObsCatch}{Array of observed catches dimensioned by n_regions, n_years, n_fish_fleets}

\item{Catch_Type}{Matrix of catch types dimensioned n_years, n_fish_fleets, == 0 indicates years to use aggregated catch, == 1 indicates years to use region specific catch}

\item{UseCatch}{Array of indicator variables dimensioned by n_regions, n_years, n_fish_fleets, == 0 don't use catch, == 1 use catch and fit}

\item{Use_F_pen}{Indicator number specifying whether a f penalty is applied == 0, don't apply, == 1 apply}

\item{est_all_regional_F}{Indicator number specifying whether all regional fishing mortality deviates are estimated, == 0 some F's are not regional and are aggregated, == 1 all fishing mortality deviates are regional}

\item{Catch_Constant}{A vector dimensioned by n_fish_fleets specifying the constant to add to catch observations}

\item{sigmaC_spec}{Observation error specificaiotn for catch. Default behavior is to fix it for all regions and fishery fleets. Other options include: "est_shared_f" which estimates it but shares it across fishery fleets, "est_shared_r" which estimates it but shares across regions (unique for each fleet), "est_shared_r_f which estimates it but shares across regions and fleets, "fix" which fixes it at the starting value, and "est_all", which estimates them all}

\item{sigmaF_spec}{Process error specificaiotn for fishing mortality Default behavior is to fix it for all regions and fishery fleets. Other options include: "est_shared_f" which estimates it but shares it across fishery fleets, "est_shared_r" which estimates it but shares across regions (unique for each fleet), "est_shared_r_f which estimates it but shares across regions and fleets, "fix" which fixes it at the starting value, and "est_all", which estimates them all}

\item{sigmaF_agg_spec}{Process error specificaiotn for fishing mortality when aggregated. Default behavior is to fix it for all fishery fleets. Other options include: "est_shared_f" which estimates it but shares it across fishery fleets, "fix" which fixes it at the starting value, and "est_all", which estimates them all}

\item{...}{Additional arguments specifying starting values for ln_sigmaC and ln_sigmaF}
}
\description{
Setup fishing mortality and catch observations
}
