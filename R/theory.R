#' @keywords internal
get_lprior <- function(mu, tau, sei) {
  e_fisher_i <- function(se) {
    si <- sqrt(tau ^ 2 + se ^ 2)

    kmm <- -si ^ (-2)
    kms <- 0
    kss <- -2 * tau ^ 2 * si ^ (-4)

    matrix(c(-kmm, -kms, -kms, -kss), nrow = 2, ncol = 2)
  }

  e_fisher <- purrr::map(sei, e_fisher_i) |> purrr::reduce(`+`)
  log(sqrt(det(e_fisher)))
}

#' @keywords internal
get_nll <- function(mu, tau, yi, sei) {
  si <- sqrt(tau ^ 2 + sei ^ 2)
  sum(log(si * sqrt(2 * pi)) + 0.5 * si ^ (-2) * (yi - mu) ^ 2)
}

#' @keywords internal
nlpost <- function(mu, tau, yi, sei) {
  joint_nll <- get_nll(mu, tau, yi, sei) # negative log-likelihood
  joint_lprior <- get_lprior(mu, tau, sei) # log-prior
  joint_nll - joint_lprior # log-posterior
}

#' @keywords internal
mle_params <- function(mu_start, tau_start, yi, sei) {
  nlpost_fun <- function(mu, tau) nlpost(mu, tau, yi, sei)
  stats4::mle(minuslogl = nlpost_fun,
              start = list(mu = mu_start, tau = tau_start),
              method = "Nelder-Mead")
}
