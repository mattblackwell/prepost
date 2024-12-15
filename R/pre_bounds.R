## PRE BOUNDS ----

#' Run pre-treatment bounds.
#'
#' @inheritParams prepost_gibbs
#' @param conf_level A numeric indicating the confidence level for the bootstrap
#'   confidence intervals.
#' @param outcome_mono A integer indicating the direction of the priming
#'   monotonicity assumption. The default value `1` indicates that asking the
#'   moderator question in the pre-test moves outcomes in a positive direction
#'   for all units. The value `-1` indicates it moves outcomes in a negative
#'   direction for all units.
#'
#' @return A list object containing bounds.
#'
#' @examples
#' data(delponte)
#' pre_bounds(
#'   formula = angry_bin ~ t_commonality,
#'    data = delponte,
#'   moderator = ~ itaid_bin
#' )
#' @importFrom stats reformulate
#' @export
pre_bounds <- function(formula, data, moderator,
                        conf_level = 0.95,
                        outcome_mono = 1L) {


  # Extract, outcome, moderator, treatment variables from formulas

  outcome = all.vars(formula)[1]
  treat = all.vars(formula)[2]
  moderator = all.vars(moderator)[1]

  form <- reformulate(c(treat, moderator), response = outcome)
  mf <- model.frame(form, data)

  Y <- as.numeric(unlist(mf[, outcome]))
  D <- as.numeric(unlist(mf[, moderator]))
  T <- as.numeric(unlist(mf[, treat]))
  Z <- rep(0, length(Y))


  pre_est <- mean(Y[T == 1 & D == 1]) - mean(Y[T == 0 & D == 1]) -
    mean(Y[T == 1 & D == 0]) + mean(Y[T == 0 & D == 0])

  out <- list()
  p <- tapply(Y, interaction(T, Z, D, sep = ""), mean)
  q <- tapply(rep(1, length(Y)), interaction(T, Z, D, sep = ""), sum)
  s2 <- p * (1 - p) / q

  if (outcome_mono == 1) {
    s_upper <- sqrt(s2["101"] + s2["000"])
    s_lower <- sqrt(s2["100"] + s2["001"])
    out$lower <- -(p["100"] + p["001"])
    out$upper <- p["101"] + p["000"]
  } else if (outcome_mono == -1) {
    s_upper <- sqrt(s2["001"] + s2["100"])
    s_lower <- sqrt(s2["101"] + s2["000"])
    out$lower <- p["101"] + p["000"] - 2
    out$upper <- 2 - p["001"] - p["100"]
  }

  names(out$lower) <- ""
  names(out$upper) <- ""

  im <- imb.man.ci(out$lower, out$upper, s_lower, s_upper,
                        N = length(Y), alpha = conf_level)
  out$ci_lower <- unname(im[1])
  out$ci_upper <- unname(im[2])

  out$pre_est <- pre_est
  return(out)
}



#' Run sensitivity analysis on pre-test design
#'
#' @inheritParams pre_bounds
#' @param t_by Numeric indicating the grid spacing for the
#' \eqn{\theta} parameter that restricts what proportion of units have
#' their outcomes affected by the pre vs post-measurement of the
#' moderator.
#'
#' @return A list object containing sensitivity output.
#'
#' @examples
#' pre_sens(formula = angry_bin ~ t_commonality,
#'   data = delponte,
#'   moderator = ~ itaid_bin,
#'   t_by = 0.1
#' )
#' @export
pre_sens <- function(formula, data,  moderator,
                      t_by = 0.05, conf_level = 0.95,
                      outcome_mono = 1L) {

  


  outcome = all.vars(formula)[1]
  treat = all.vars(formula)[2]
  moderator = all.vars(moderator)[1]

  form <- reformulate(c(treat, moderator), response = outcome)
  mf <- model.frame(form, data)

  Y <- as.numeric(unlist(mf[, outcome]))
  D <- as.numeric(unlist(mf[, moderator]))
  T <- as.numeric(unlist(mf[, treat]))
  Z <- rep(0, length(Y))


  p <- tapply(Y, interaction(T, Z, D, sep = ""), mean)
  q <- tapply(rep(1, length(Y)), interaction(T, Z, D, sep = ""), sum)
  s2 <- p * (1 - p) / q  
  
  pre_est <- p["101"] - p["001"] - p["100"] + p["000"]  
  pre_se <- sqrt(s2["101"] + s2["001"] + s2["100"] + s2["000"])
  
  if (outcome_mono == 1) {
    t_max <- max(p["101"], p["100"], p["001"], p["000"])
  } else {
    t_max <- 1 - min(p["101"], p["100"], p["001"], p["000"])
  }
  thetas <- seq(0, t_max, by = t_by)

  out <- list()
  out$thetas <- thetas
  out$lower <- pre_est - 2 * thetas
  out$upper <- pre_est + 2 * thetas

  
  out$ci_lower <- rep(NA, length(thetas))
  out$ci_upper <- rep(NA, length(thetas))

  for (k in seq_along(thetas)) {
    im <- imb.man.ci(out$lower[k], out$upper[k], pre_se, pre_se,
                     N = length(Y), alpha = conf_level)
    out$ci_lower[k] <- unname(im[1])
    out$ci_upper[k] <- unname(im[2])

  }
  return(out)
  

}
