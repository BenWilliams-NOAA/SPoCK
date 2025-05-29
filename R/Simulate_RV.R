#' Simulate logistic normal variables
#'
#' @param exp Expected proportions (should sum to 1)
#' @param pars Parameters for a logistic normal:
#'   - comp_like = 2: pars = c(sigma)
#'   - comp_like = 3: pars = c(sigma, rho_age)
#'   - comp_like = 4: pars = c(sigma, rho_age, rho_sex)
#'   - comp_like = 5: pars = c(sigma, rho_age, rho_sex, rho_region)
#'    (iid = 1 parameter; AR1 = 2 parameters; 2D, by age and sex = 3 parameters; 3D, by age, sex, and region = 4 parameters)
#' @param comp_like Likelihood structure (2 = iid, 3 = ar1, 4 = 2d, 5 = 3d)
#' @importFrom MASS mvrnorm
#' @keywords internal
#'
rlogistnormal <- function(exp, pars, comp_like) {

  # Remove last bin and compute log-ratios
  mu <- log(exp[-length(exp)]) - log(exp[length(exp)])

  # Determine covariance structure
  Sigma <- switch(comp_like,
                  `2` = diag(pars[1]^2, length(mu)),

                  `3` = {
                    Sigma_raw <- get_AR1_CorrMat(n_ages, pars[2]) * pars[1]^2
                    Sigma_raw[-nrow(Sigma_raw), -ncol(Sigma_raw)]
                  },

                  `4` = {
                    Sigma_raw <- kronecker(
                      get_AR1_CorrMat(n_ages, pars[2]),
                      get_Constant_CorrMat(n_sexes, pars[3])
                      ) * pars[1]^2
                    Sigma_raw[-nrow(Sigma_raw), -ncol(Sigma_raw)]
                  },

                  `5` = {
                    Sigma_raw <- kronecker(
                      kronecker(
                        get_AR1_CorrMat(n_ages, pars[2]),
                        get_Constant_CorrMat(n_sexes, pars[3])
                      ),
                      get_Constant_CorrMat(n_regions, pars[4])
                    ) * pars[1]^2
                    Sigma_raw[-nrow(Sigma_raw), -ncol(Sigma_raw)]
                  },

                  stop("Unsupported comp_like type")
  )

  # simulate and transform back to proportions
  x <- MASS::mvrnorm(1, mu, Sigma)
  p <- exp(x)/(1 + sum(exp(x)))
  p <- c(p, 1 - sum(p))

  return(p)
}


#' Simulate dirichlet multinomial draws
#'
#' @param n Number of sims
#' @param N Sum of observations
#' @param alpha Concentration parameter
#' @importFrom stats rgamma rmultinom
#' @keywords internal
rdirM <- function(n, N, alpha) {

  # Get dirichlet draws
  rdirichlet <- function(alpha) {
    x <- stats::rgamma(length(alpha), shape = alpha, scale = 1)
    return(x / sum(x))
  }

  # Generate DM samples
  result <- replicate(n, {
    p <- rdirichlet(alpha)
    counts <- stats::rmultinom(1, size = N, prob = p)
    as.vector(counts)
  })

  return(result)
}

#' Generate Recruitment Values based on Inverse Gaussian Distribution
#'
#' @param sims Number of Simulations
#' @param recruitment Recruitment vector
#' @importFrom stats rnorm runif
#' @returns Random variables following an inverse gaussian
#' @keywords internal
rinvgauss_rec <- function(sims,
                          recruitment
                          ) {

  # calculate arithmetic and harmonic means of recruitment
  a_mean = mean(recruitment)
  h_mean = 1 / mean(1 / recruitment)

  # calculate dispersion parameter based on ratio of means
  gamma = a_mean / h_mean
  delta = 1 / (gamma - 1)
  cv = sqrt(1 / delta)

  # generate standard normal squared for transformation
  psi = stats::rnorm(sims)^2

  # inverse Gaussian transformation components
  term = sqrt(4 * delta * psi + psi^2)
  omega = a_mean * (1 + (psi - term) / (2 * delta))
  zeta  = a_mean * (1 + (psi + term) / (2 * delta))

  # probability mixture for inverse Gaussian
  gtheta = a_mean / (a_mean + omega)
  rv = ifelse(runif(sims) <= gtheta, omega, zeta)

  return(rv)

}

