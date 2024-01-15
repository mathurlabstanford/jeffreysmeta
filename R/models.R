#' Right-truncated meta-analysis
#'
#' Fits right-truncated meta-analysis (RTMA), a bias correction for the joint
#' effects of p-hacking (i.e., manipulation of results within studies to obtain
#' significant, positive estimates) and traditional publication bias (i.e., the
#' selective publication of studies with significant, positive results) in
#' meta-analyses. This method analyzes only nonaffirmative studies (i.e., those
#' with significant, positive estimates). You can pass all studies in the meta-analysis
#' or only the nonaffirmative ones; if the former, the function will still analyze only
#' the nonaffirmative ones.
#'
#' @inheritParams metabias::params
#' @param stan_control List passed to [rstan::sampling()] as the `control`
#'   argument.
#' @param parallelize Logical indicating whether to parallelize sampling.
#'
#' @return An object of class [metabias::metabias()], a list containing:
#' \describe{
#'   \item{data}{A tibble with one row per study and the columns
#'               `r meta_names_str("data")`.}
#'   \item{values}{A list with the elements `r meta_names_str("values")`.
#'                 `optim_converged` indicates whether the optimization to find
#'                 the posterior mode converged.}
#'   \item{stats}{A tibble with two rows and the columns
#'                `r meta_names_str("stats")`. We recommend reporting the `mode`
#'                for the point estimate; `median` and `mean` represent
#'                posterior medians and means respectively.}
#'   \item{fit}{A `stanfit` object (the result of fitting the RTMA model).}
#' }
#'
#' @export
#'
#' @references
#' \insertRef{mathur2022phacking}{metabias}
#'
#' @examples
#' \donttest{
#' money_priming <- jeffreys_meta(money_priming_meta$yi, money_priming_meta$vi,
#'                                     parallelize = FALSE)
#' }
jeffreys_meta <- function(yi, # data
                          vi,
                          sei,

                          ci_level = 0.95,
                          stan_control = list(adapt_delta = 0.98,
                                              max_treedepth = 20),
                          parallelize = TRUE) {

  if (missing(vi) && missing(sei)) stop("Must specify 'vi' or 'sei' argument.")
  if (missing(vi)) vi <- sei ^ 2
  if (missing(sei)) sei <- sqrt(vi)
  if (length(sei) != length(yi)) stop(
    "Length of 'vi' or 'sei' must match that of 'yi'."
  )
  if (any(sei < 0)) stop("vi or sei should never be negative.")

  k <- length(yi)

  dat <- tibble(yi = yi, vi = vi, sei = sei)
 
  stan_data <- list(y = array(yi), sei = array(sei),
                    k = k)

  vals <- list(ci_level = ci_level,
               k = k)

  if (parallelize) options(mc.cores = parallel::detectCores())
  stan_fit <- rstan::sampling(stanmodels$jeffreys_meta,
                              data = stan_data,
                              control = stan_control,
                              init = function() list(mu = 0, tau = 1))

  stan_extract <- rstan::extract(stan_fit)
  medians <- c(median(stan_extract$mu), median(stan_extract$tau))

  index_maxlp <- which.max(stan_extract$log_post)
  mu_maxlp <- stan_extract$mu[index_maxlp]
  tau_maxlp <- stan_extract$tau[index_maxlp]
  mle_fit <- mle_params(mu_maxlp, tau_maxlp, yi, sei)
  modes <- c(mle_fit@coef[["mu"]], mle_fit@coef[["tau"]])
  optim_converged <- mle_fit@details$convergence == 0
  vals$optim_converged <- optim_converged

  stan_summary <- rstan::summary(stan_fit)$summary |>
    as_tibble(rownames = "param")

  stan_ci <- function(param, q) as.numeric(quantile(stan_extract[[param]], q))
  alpha <- 1 - ci_level
  cil <- c(stan_ci("mu", alpha / 2), stan_ci("tau", alpha / 2))
  ciu <- c(stan_ci("mu", 1 - alpha / 2), stan_ci("tau", 1 - alpha / 2))

  stan_stats <- stan_summary |>
    filter(.data$param %in% c("mu", "tau")) |>
    select(.data$param, .data$mean, se = .data$se_mean,
           .data$n_eff, r_hat = .data$Rhat) |>
    mutate(mode = modes, median = medians, .after = .data$param) |>
    mutate(ci_lower = cil, ci_upper = ciu, .after = .data$se)

  metabias::metabias(data = dat, values = vals,
                     stats = stan_stats, fits = stan_fit)

}
