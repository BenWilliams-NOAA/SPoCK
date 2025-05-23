% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/do_jitter.R
\name{do_jitter}
\alias{do_jitter}
\title{Run Jitter Analysis}
\usage{
do_jitter(
  data,
  parameters,
  mapping,
  random = NULL,
  sd,
  n_jitter,
  n_newton_loops,
  do_par,
  n_cores
)
}
\arguments{
\item{data}{Data list to make obj}

\item{parameters}{Parameter list to make obj}

\item{mapping}{Mapping list to make obj}

\item{random}{Character of random effects}

\item{sd}{sd for jitter (additive)}

\item{n_jitter}{Number of jitters to do}

\item{n_newton_loops}{Number of newton loops to do}

\item{do_par}{Whether to do paralleizaiton or not (boolean)}

\item{n_cores}{Number of cores to use}
}
\value{
Dataframe of jitter values
}
\description{
Run Jitter Analysis
}
\examples{
\dontrun{
   library(ggplot2)
   # get jitter values
   jit <- do_jitter(data = data,
                 parameters = parameters,
                 mapping = mapping,
                 random = NULL,
                 sd = 0.1,
                 n_jitter = 100,
                 n_newton_loops = 3,
                 do_par = TRUE,
                 n_cores = 8)

   # get proportion converged
   prop_converged <- jit \%>\%
   filter(Year == 1, Type == 'Recruitment') \%>\%
     summarize(prop_conv = sum(Hessian) / length(Hessian))

   # get final model results
   final_mod <- reshape2::melt(sabie_rtmb_model$rep$SSB) \%>\% rename(Region = Var1, Year = Var2) \%>\%
   mutate(Type = 'SSB') \%>\%
     bind_rows(reshape2::melt(sabie_rtmb_model$rep$Rec) \%>\%
     rename(Region = Var1, Year = Var2) \%>\% mutate(Type = 'Recruitment'))

   # comparison of SSB and recruitment
  ggplot() +
    geom_line(jit, mapping = aes(x = Year + 1959, y = value, group = jitter, color = Hessian), lwd = 1) +
    geom_line(final_mod, mapping = aes(x = Year + 1959, y = value), color = "black", lwd = 1.3 , lty = 2) +
    facet_grid(Type~Region, scales = 'free',
               labeller = labeller(Region = function(x) paste0("Region ", x),
                                   Type = c("Recruitment" = "Age 2 Recruitment (millions)", "SSB" = 'SSB (kt)'))) +
    labs(x = "Year", y = "Value") +
    theme_bw(base_size = 20) +
    scale_color_manual(values = c("red", 'grey')) +
    geom_text(data = jit \%>\% filter(Type == 'SSB', Year == 1, jitter == 1),
              aes(x = Inf, y = Inf, label = paste("Proportion Converged: ", round(prop_converged$prop_conv, 3))),
              hjust = 1.1, vjust = 1.9, size = 6, color = "black")

   # compare jitter of max gradient and hessian PD
   ggplot(jit, aes(x = jitter, y = jnLL, color = Max_Gradient, shape = Hessian)) +
     geom_point(size = 5, alpha = 0.3) +
     geom_hline(yintercept = min(sabie_rtmb_model$rep$jnLL), lty = 2, size = 2, color = "blue") +
     facet_wrap(~Hessian, labeller = labeller(
       Hessian = c("FALSE" = "non-PD Hessian", "TRUE" = 'PD Hessian')
     )) +
     scale_color_viridis_c() +
     theme_bw(base_size = 20) +
     theme(legend.position = "bottom") +
     guides(color = guide_colorbar(barwidth = 15, barheight = 0.5)) +
     labs(x = 'Jitter') +
     geom_text(data = jit \%>\% filter(Hessian == TRUE, Year == 1, jitter == 1),
               aes(x = Inf, y = Inf, label = paste("Proportion Converged: ", round(prop_converged$prop_conv, 3))),
               hjust = 1.1, vjust = 1.9, size = 6, color = "black")
}
}
