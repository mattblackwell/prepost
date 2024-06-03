



imb.man.ci <- function(lo, hi, lo.sd, hi.sd, N, alpha = 0.95) {
  delta <- sqrt(N) * (hi - lo) / (max(lo.sd, hi.sd))
  Cs <- seq(0, 10, by = 0.001)
  Cseq <- abs(pnorm(Cs + delta) - pnorm(-Cs) - alpha)
  Cn <- Cs[Cseq == min(Cseq)]

  ci.lo <- lo - Cn * lo.sd
  ci.hi <- hi + Cn * hi.sd
  return(ci = c(ci.lo, ci.hi))
}


compute_strata_probs <- function(Y, D, T, Z) {


  P_levs <- list(Y1 = c(1, 0), D1 = c(1, 0), T = c(1, 0))

  ## create names of P strata to be able to grep
  ## P_TZD
  P_vals <- expand.grid(P_levs)
  nms <- do.call(paste0, P_vals)
  nstrata <- length(nms)
  P1 <- rep(0, times = nstrata)
  names(P1) <- nms
  Ptab <- table(paste0(Y, D, T)[Z == 1])
  P1[names(Ptab)] <- Ptab
  P_den <- ave(P1, P_vals[, c("T"), drop = FALSE], FUN = sum)
  if (sum(P1) != 0) {
    P <- P1 / P_den
  } else {
    P <- P1
  }

  ## V_YTD
  V_levs <- list(Y0 = c(1, 0), T = c(1, 0), D0 = c(1, 0))
  V_vals <- expand.grid(V_levs)
  V_nms <- do.call(paste0, V_vals)
  V_nstrata <- length(V_nms)
  V1 <- rep(0, times = V_nstrata)
  names(V1) <- V_nms
  Vtab <- table(paste0(Y, T, D)[Z == 0])
  V1[names(Vtab)] <- Vtab
  V_den <- ave(V1, V_vals[, c("T", "D0")], FUN = sum)
  if (sum(V1) != 0) {
    V <- V1 / V_den
  } else {
    V <- V1
  }


  if (sum(Z == 0) > 0) {
    Q <- mean(D[Z == 0])
  } else {
    ## if we only have post, calculate Q under stable mod under
    ## control
   Q <- mean(D[Z == 1 & T == 0])
  }

  ## calculate empirical direction of the placement effect on the moderator
  MM <- rep(NA, times = 2)
  OM <- rep(NA, times = 2)
  if (sum(Z == 1) > 0) {
    for (j in c(0, 1)) {
      dd_z1 <- mean(D[Z == 1 & T == j])
      dd_z0 <- mean(D[Z == 0])
      MM[j + 1] <- sign(dd_z1 - dd_z0)
      names(MM)[j + 1] <- j

      yy_z1 <- mean(Y[Z == 1 & T == j])
      yy_z0 <- mean(Y[Z == 0])
      OM[j + 1] <- sign(yy_z1 - yy_z0)
      names(OM)[j + 1] <- j
    }

  }



  return(list(P = P, P_den = P_den, V = V, V_den = V_den,
              Q = Q, MM = MM, OM = OM))
}



## POST BOUNDS ----

#' Run post-treatment bounds.
#'
#' @inheritParams prepost_gibbs
#' @param sims An integer indicating the number of simulations for the
#' bootstrap confidence intervals for the bounds.
#' @param conf_level A numeric indicating the confidence level for the
#' bootstrap confidence intervals.
#' @param moderator_mono A integer or vector of length 2 indicating
#' if the bounds should assume monotonicity of the effect of the
#' post-test on the moderator with `1` indicating that the post-test
#' effect is positive and `-1` indicating that it is negative. The
#' vector of length 2 allows the monotonicity assumption to vary by
#' treatment status with the first entry being for control and the
#' second for treated.
#' @param stable_mod A logical value indicating if the bounds should
#' assume that the moderator is unaffected by pre-vs-post measurement
#' under the control condition.
#' @param nondiff A logical value indicating if the bounds should
#' assume the treatment effect on the moderator is independent of the
#' potential outcomes.
#' @param progress A logical indicating if progress bars should be
#' displayed. Defaults to TRUE.
#'
#' @return A list object containing bounds.
#'
#' @examples
#' data(delponte)
#' post_bounds(
#'   formula = angry_bin ~ t_commonality,
#'    data = delponte,
#'   moderator = ~ itaid_bin
#' )
#' @importFrom stats pnorm sd optimize
#' @export
post_bounds <- function(formula, data, moderator, sims = 1000,
                        conf_level = 0.95,
                        moderator_mono = NULL, stable_mod = FALSE,
                        nondiff = FALSE, progress = TRUE) {


  # Extract, outcome, moderator, treatment variables from formulas

  outcome = all.vars(formula)[1]
  treat = all.vars(formula)[2]
  moderator = all.vars(moderator)[1]

  # SH: coercing these vars to the correct format due to issue with tibbles
  Y <- as.numeric(unlist(data[, outcome]))
  D <- as.numeric(unlist(data[, moderator]))
  T <- as.numeric(unlist(data[, treat]))
  Z <- rep(1, length(Y))


  post_est <- mean(Y[T == 1 & D == 1]) - mean(Y[T == 0 & D == 1]) -
    mean(Y[T == 1 & D == 0]) + mean(Y[T == 0 & D == 0])


  q <- tapply(D, interaction(T, Z, sep = ""), mean)
  if (is.logical(moderator_mono)) {
    if (moderator_mono == TRUE & stable_mod) {
      if (q["11"] > q["01"]) {
        moderator_mono <- c(1, 1)
      } else {
        moderator_mono <- c(-1,-1)
      }
    }
  }

  if (length(moderator_mono) > 0) {
    if (length(moderator_mono) == 1) {
      moderator_mono <- rep(moderator_mono, 2)
    }

    if (diff(moderator_mono) != 0) {
      if (sign(diff(rev(moderator_mono))) != sign(q["11"] - q["01"])) {
        stop("moderator monotonicity assumptions not compatible with data",
             call. = FALSE)
      }
    }
  }

  if (identical(moderator_mono, c(1, 1))) {
    if (stable_mod & q["11"] < q["01"]) {
      stop("stable moderator and monotonicity not compatible with data",
           call. = FALSE)
    }
  }
  if (identical(moderator_mono, c(-1, -1))) {
    if (stable_mod & q["11"] > q["01"]) {
      stop("stable moderator and monotonicity not compatible with data",
           call. = FALSE)
    }
  }

  out <- list()
  p <- tapply(Y, interaction(T, Z, D, sep = ""), mean)
  q <- tapply(D, interaction(T, Z, sep = ""), mean)
  bounds <- post_factory(p, q, moderator_mono, stable_mod)
  out$lower <- unname(bounds$L)
  out$upper <- unname(bounds$U)

  if (nondiff) {
    nd_out <- post_bounds_nondiff(p, q)
    out$lower <- nd_out$lower
    out$upper <- nd_out$upper
  }

  boot_lo <- rep(NA, times = sims)
  boot_hi <- rep(NA, times = sims)
  if (progress) cat("Bootstrap running...\n")
  for (b in 1:sims) {
    if ((100 * b / sims) %% 10 == 0 & progress) {
      pb <-  c(rep("==", times = (10 * b / sims)),
               rep("  ", times = 10 - (10 * b / sims)))
      pb <- paste0(c("0% |", pb, "| 100%\n"), collapse = "")

      cat(pb)
    }
    star <- sample(seq_len(nrow(data)), size = nrow(data), replace = TRUE)
    Ys <- Y[star]
    Ds <- D[star]
    Ts <- T[star]
    Zs <- Z[star]

    p <- tapply(Ys, interaction(Ts, Zs, Ds, sep = ""), mean)
    q <- tapply(Ds, interaction(Ts, Zs, sep = ""), mean)
    bounds <- post_factory(p, q, moderator_mono, stable_mod)
    boot_hi[b] <- bounds$L
    boot_lo[b] <- bounds$U

    if (nondiff) {
      nd_out <- post_bounds_nondiff(p, q)
      boot_hi[b] <- nd_out$lower
      boot_lo[b] <- nd_out$upper
    }
  }

  boot_sd_lo <- sd(boot_lo, na.rm = TRUE)
  boot_sd_hi <- sd(boot_hi, na.rm = TRUE)
  boot_ci <- imb.man.ci(out$lower, out$upper, boot_sd_lo, boot_sd_hi,
                        N = nrow(data), alpha = conf_level)
  out$ci_lower <- unname(boot_ci[1])
  out$ci_upper <- unname(boot_ci[2])

  out$post_est <- post_est
  return(out)
}

post_factory <- function(p, q, moderator_mono, stable_mod) {
  core <- function(Q0) {
    core <- p["111"] * (q["11"] / Q0) - p["011"] * (q["01"] / Q0)
    core <- core - p["110"] * (1 - q["11"]) / (1 - Q0)
    core <- core + p["010"] * (1 - q["01"]) / (1 - Q0)
  }
  if (!is.null(moderator_mono)) {
    if (identical(moderator_mono, c(1, 1))) {
      upper <- function(Q0) {
        den <- (Q0 * (1 - Q0))
        core(Q0) + pmin(Q0 - p["111"] * q["11"], 0) / den +
          pmin(p["011"] * q["01"], q["01"] - Q0) / den
      }
      U <- max(upper(c(
        p["111"] * q["11"], q["01"]  * (1 - p["011"]),
        min(q["11"], q["01"])
      )))
      lower <- function(Q0) {
        den <- (Q0 * (1 - Q0))
        core(Q0) - pmin(p["111"] * q["11"], q["11"] - Q0) / den +
          pmax((p["011"] * q["01"] - Q0), 0) / den
      }
      L <- min(lower(c(
        p["011"] * q["01"], q["11"] * (1 - p["111"]),
        min(q["11"], q["01"])
      )))
    } else if (identical(moderator_mono, c(1, -1))) {
      upper <- function(Q0) {
        den <- Q0 * (1 - Q0)
        core(Q0) - pmax(p["111"] * q["11"] - Q0, 0) / den -
          pmax(p["010"] * (1 - q["01"]) - (1 - Q0), 0) / den
      }
      U <- max(upper(c(
        p["111"] * q["11"], 1 - p["010"] * (1 - q["01"]),
        q["01"], q["11"]
      )))
      lower <- function(Q0) {
        den <- Q0 * (1 - Q0)
        core(Q0) - pmin(p["111"] * q["11"], q["11"] - Q0) / den -
          pmin(p["010"] * (1 - q["01"]), Q0 - q["01"]) / den
      }
      L <- min(lower(c(
        q["11"] * (1 - p["111"]),
        p["010"] * (1 - q["01"]) + q["01"],
        q["01"], q["11"]
      )))
    } else if (identical(moderator_mono, c(-1, 1))) {
      upper <- function(Q0){
        den <- Q0 * (1 - Q0)
        core(Q0) + pmin(p["011"] * q["01"], q["01"] - Q0) / den +
           pmin(p["110"] * (1 - q["11"]), Q0 - q["11"]) / den
      }
      U <- max(upper(c(
        q["01"] * (1 - p["011"]),
        p["110"] * (1 - q["11"]) + q["11"],
        q["01"], q["11"]
      )))
      lower <- function(Q0) {
        den <- Q0 * (1 - Q0)
        core(Q0) + pmax(p["011"] * q["01"] - Q0, 0) / den +
          pmax(p["110"] * (1 - q["11"]) - (1 - Q0), 0) / den

      }
      L <- min(lower(c(
        p["011"] * q["01"],
        1 - p["110"] * (1 - q["11"]),
        q["01"], q["11"]
      )))
    } else if (identical(moderator_mono, c(-1, -1))) {
      upper <- function(Q0) {
        den <- Q0 * (1 - Q0)
        core(Q0) + pmin(p["110"] * (1 - q["11"]), Q0 - q["11"]) / den -
          pmax(p["010"] * (1 - q["01"]) - (1 - Q0), 0) / den
      }
      U <- max(upper(c(
        q["11"] + p["110"] * (1 - q["11"]),
        1 - p["010"] * (1 - q["01"]),
        max(q["11"], q["01"])
      )))
      lower <- function(Q0) {
        den <- Q0 * (1 - Q0)
        core(Q0) + pmax(p["110"] * (1 - q["11"]) - (1 - Q0), 0) / den -
          pmin(p["010"] * (1 - q["01"]), Q0 - q["01"]) / den
      }
      L <- min(lower(c(
        1 - p["110"] * (1 - q["11"]),
        q["01"] + p["010"] * (1 - q["01"]),
        max(q["11"], q["01"])
      )))
    }
    if (stable_mod) {
      U <- upper(q["01"])
      L <- lower(q["01"])
    }
  } else if (stable_mod) {
    den <- (q["01"] * (1 - q["01"]))
    U <- core(q["01"]) +
      pmin(
        p["110"] * (1 - q["11"]) / den,
        1 / (1 - q["01"]) - p["111"] * q["11"] / den
      )
    L <- core(q["01"]) +
      pmax(
        -p["111"] * q["11"] / den,
        (p["110"] * (1 - q["11"]) - (1 - q["01"])) / den
      )
  } else {
    P1 <- p["111"] * q["11"] + p["110"] * (1 - q["11"])
    P0 <- p["011"] * q["01"] + p["010"] * (1 - q["01"])
    upper <- function(Q0) {
      pmin(
        2,
        (1 + P1 - P0) / Q0,
        (1 - P1 + P0) / (1 - Q0),
        (1 - P1) / (1 - Q0) + (1 - P0) / Q0,
        (P1) / Q0 + (P0) / (1 - Q0)
      )
    }

    lower <- function(Q0) {
      pmax(
        -2,
        -(1 - P1 + P0) / Q0,
        -(1 + P1 - P0) / (1 - Q0),
        -(P1) / (1-Q0) - (P0) / Q0,
        -(1 - P1) / Q0 - (1 - P0) / (1 - Q0)
      )
    }
    U <- max(upper(c(P1, 1 - P0)))
    L <- min(lower(c(1 - P1, P0)))
  }

  return(list(U = U, L = L))


}

## Sensitivity analysis ----

#' Run sensitivity analysis on post-measurement design
#'
#' @inheritParams post_bounds
#' @param g_by Numeric indicating the grid spacing for the
#' \eqn{\gamma} parameter that places an upper bound on the proportion
#' of units whose moderator is affected by treatment.
#' @param q_by Numeric indicating the grid spacing for the mean of the
#' moderator under a pre-test measurement.
#' @param g_max Numeric indicating the maximum value of the \eqn{\gamma} parameter. 
#'
#' @return A list object containing sensitivity output.
#'
#' @examples
#' post_sens(formula = angry_bin ~ t_commonality,
#'   data = delponte,
#'   moderator = ~ itaid_bin,
#'   g_by = 0.1
#' )
#' @export
post_sens <- function(formula, data,  moderator,
                      g_by, g_max = 1, q_by, sims = 1000, conf_level = 0.95,
                      moderator_mono = NULL, stable_mod = FALSE,
                      progress = TRUE) {


  # Extract, outcome, moderator, treatment, and covariates variables from formulas

  outcome = all.vars(formula)[1]
  treat = all.vars(formula)[2]
  moderator = all.vars(moderator)[1]


  Y <- data[, outcome]
  D <- data[, moderator]
  T <- data[, treat]
  Z <- rep(1, length(Y))

  obs <- compute_strata_probs(Y = Y, D = D, T = T, Z = Z)
  p <- tapply(Y, interaction(T, Z, D, sep = ""), mean)
  q <- tapply(D, interaction(T, Z, sep = ""), mean)

  min_gamma <- abs(q["11"] - q["01"])
  gammas <- seq(min_gamma, g_max, by = g_by)

  sens_out <- list()
  sens_out$gamma <- gammas
  sens_out$lower <- rep(NA, times = length(gammas))
  sens_out$upper <- rep(NA, times = length(gammas))

  for (g in seq_along(gammas)) {
    q_min <- max(0, q["01"] - gammas[g], q["11"] - gammas[g])
    q_max <- min(1, gammas[g] + q["01"], gammas[g] + q["11"])
    if (missing(q_by)) {
      upper <- optimize(gamma_bounds, interval = c(q_min, q_max), maximum = TRUE,
                        p = p, q = q, gamma = gammas[g], hi = TRUE)
      lower <- optimize(gamma_bounds, interval = c(q_min, q_max), maximum = FALSE,
                        p = p, q = q, gamma = gammas[g], hi = FALSE)
      sens_out$upper[g] <- upper$objective
      sens_out$lower[g] <- lower$objective
    } else {
      Q0 <- seq(q_min, q_max, by = q_by)
      uppers <- gamma_bounds(Q0, p = p, q = q, gamma = gammas[g], hi = TRUE)
      lowers <- gamma_bounds(Q0, p = p, q = q, gamma = gammas[g], hi = FALSE)
      sens_out$upper[g] <- max(uppers)
      sens_out$lower[g] <- min(lowers)
    }
  }


  sens_lo <- sens_hi <- array(NA, c(sims, length(sens_out$gamma)))
  for (b in 1:sims) {
    if ((100 * b / sims) %% 10 == 0 & progress) {
      pb <-  c(rep("==", times = (10 * b / sims)),
               rep("  ", times = 10 - (10 * b / sims)))
      pb <- paste0(c("0% |", pb, "| 100%\n"), collapse = "")

      cat(pb)
    }
    star <- sample(seq_len(nrow(data)), size = nrow(data), replace = TRUE)
    Ys <- Y[star]
    Ds <- D[star]
    Ts <- T[star]
    Zs <- Z[star]
    ostar <- compute_strata_probs(Y = Ys, D = Ds, T = Ts, Z = Zs)

    ps <- tapply(Ys, interaction(Ts, Zs, Ds, sep = ""), mean)
    qs <- tapply(Ds, interaction(Ts, Zs, sep = ""), mean)

    for (g in seq_along(gammas)) {
      q_min <- max(0, qs["01"] - gammas[g], qs["11"] - gammas[g])
      q_max <- min(1, gammas[g] + qs["01"], gammas[g] + qs["11"])
      if (q_min < q_max) {
        if (missing(q_by)) {
          upper <- optimize(gamma_bounds, interval = c(q_min, q_max), maximum = TRUE,
                            p = ps, q = qs, gamma = gammas[g], hi = TRUE)
          lower <- optimize(gamma_bounds, interval = c(q_min, q_max), maximum = FALSE,
                            p = ps, q = qs, gamma = gammas[g], hi = FALSE)
          sens_hi[b, g] <- upper$objective
          sens_lo[b, g] <- lower$objective
        } else {
          Q0 <- seq(q_min, q_max, by = q_by)
          uppers <- gamma_bounds(Q0, p = ps, q = qs, gamma = gammas[g], hi = TRUE)
          lowers <- gamma_bounds(Q0, p = ps, q = qs, gamma = gammas[g], hi = FALSE)
          sens_hi[b, g] <- max(uppers)
          sens_lo[b, g] <- min(lowers)
        }
      }
    }

  }
  sd_lo <- apply(sens_lo, 2, sd, na.rm = TRUE)
  sd_hi <- apply(sens_hi, 2, sd, na.rm = TRUE)
  sens_ci_lo <- rep(NA, times = length(sens_out$gamma))
  sens_ci_hi <- rep(NA, times = length(sens_out$gamma))
  for (i in seq_along(gammas)) {
    if (!is.na(sens_out$lower[i])) {
      imbman <- imb.man.ci(sens_out$lower[i], sens_out$upper[i],
                           sd_lo[i], sd_hi[i], N = nrow(data),
                           alpha = conf_level)
      sens_ci_lo[i] <- imbman[1]
      sens_ci_hi[i] <- imbman[2]
    }
  }

  sens_out$ci_lo <- sens_ci_lo
  sens_out$ci_hi <- sens_ci_hi

  return(sens_out)

}


gamma_bounds <- function(Q0, p, q, gamma, hi = TRUE) {
  rho_011_001_max <- pmin(
    1 - q["11"], Q0,
    pmax(gamma - q["01"] + Q0, 0),
    pmax(0.5 * (gamma - q["11"] + Q0), 0),
    pmax(gamma + q["01"] - q["11"], 0)
  )
  rho_110_010_max <- pmin(
    1 - Q0, q["01"],
    pmax(gamma + q["11"] - Q0, 0),
    pmax(0.5 * (gamma + q["01"] - Q0), 0),
    pmax(gamma - q["11"] + q["01"], 0)
  )
  rho_001_101_max <- pmin(
    1 - q["01"], Q0,
    pmax(gamma - q["11"] + Q0, 0),
    pmax(0.5 * (gamma - q["01"] + Q0), 0),
    pmax(gamma - q["01"] + q["11"], 0)
  )
  rho_110_100_max <- pmin(
    1 - Q0, q["11"],
    pmax(gamma + q["01"] - Q0, 0),
    pmax(0.5 * (gamma + q["11"] - Q0), 0),
    pmax(gamma - q["01"] + q["11"], 0)
  )

  rho_upper_max <- pmin(
    gamma - q["11"] + q["01"], gamma,
    q["01"] + Q0, 2 - Q0 - q["11"]
  )
  rho_lower_max <- pmin(
    gamma - q["01"] + q["11"], gamma,
    q["11"] + Q0, 2 - Q0 - q["01"]
  )

  core <- p["111"] * (q["11"] / Q0) -
    p["011"] * (q["01"] / Q0) -
    p["110"] * ((1 - q["11"]) / ( 1 - Q0)) +
    p["010"] * ((1 - q["01"]) / (1 - Q0))

  den <-  Q0 * (1 - Q0)

  if (hi) {
    bd <- core +
      pmin(
        1 / (1 - Q0) - p["111"] * q["11"] / den + p["011"] * q["01"] / den,
        1 / (1 - Q0) - p["111"] * q["11"] / den + rho_110_010_max / den,
        1 / (1 - Q0) - p["111"] * q["11"] / den + 1 / Q0 - p["010"] * (1 - q["01"]) / den,
        rho_011_001_max / den + p["011"] * q["01"] / den,
        rho_upper_max / den,
        rho_011_001_max / den + 1 / Q0 - p["010"] * (1 - q["01"]) / den,
        p["110"] * (1 - q["11"]) / den +  p["011"] * q["01"] / den,
        p["110"] * (1 - q["11"]) / den + rho_110_010_max / den,
        p["110"] * (1 - q["11"]) / den + 1 / Q0 - p["010"] * (1 - q["01"]) / den
      )
  } else {
    bd <- core +
      pmax(
        -p["111"] * q["11"] / den - p["010"] * (1 - q["01"]) / den,
        -p["111"] * q["11"] / den - rho_001_101_max / den,
        -p["111"] * q["11"] / den - 1 / (1 - Q0) + p["011"] * q["01"] / den,
        -rho_110_100_max / den - p["010"] * (1 - q["01"]) / den,
        -rho_lower_max / den,
        -rho_110_100_max / den - 1 / (1 - Q0) + p["011"] * q["01"] / den,
        -1/Q0 + p["110"] * (1 - q["11"]) / den - p["010"] * (1 - q["01"]) / den,
        -1/Q0 + p["110"] * (1 - q["11"]) / den - rho_001_101_max / den,
        -1/Q0 + p["110"] * (1 - q["11"]) / den - 1 / (1 - Q0) + p["011"] * q["01"] / den
      )
    }
  return(bd)
}


post_bounds_nondiff <- function(p, q, q_by) {
  d1 <- q["11"]
  d0 <- q["01"]
  Qmax <- min(q["01"], q["11"])
  Qmin <- max(
  (p["011"] - q["01"] * p["010"] / (1 - q["01"])) / (1 - p["010"] / (1 - q["01"])),
  (p["111"] - q["11"] * p["110"] / (1 - q["11"])) / (1 - p["110"] / (1 - q["11"]))
  )

  nd_est <- function(Q0) {
    (p["111"] + p["110"] - p["110"] * (1 - Q0) / (1 - q["11"])) / Q0 -
      (p["011"] - p["010"] * (q["01"] - Q0) / (1 - q["01"])) / Q0 -
      p["110"] / (1 - q["11"]) +
      p["010"] / (1 - q["01"])
  }

  q_upper <- nd_est(Qmax)
  q_lower <- nd_est(Qmin)

  upper <- max(q_upper, q_lower)
  lower <- min(q_upper, q_lower)
  return(list(upper = upper, lower = lower))
}

## PREPOST BOUNDS ----

#' Run Prepost bounds
#'
#' @inheritParams post_bounds
#' @param prepost A one-sided formula with syntax ~ z, where z is the indicator variable for whether the moderator was measured pre- or post-treatment. 
#' @param sims An integer specifying the number of simulations for the sensitivity analysis.
#' @param conf_level A numeric specifying level for the confidence intervals.
#' @param outcome_mono A integer or vector of length 2 indicating
#' if the bounds should assume monotonicity of the effect of the
#' post-test on the outcome with `1` indicating that the post-test
#' effect is positive and `-1` indicating that it is negative. The
#' vector of length 2 allows the monotonicity assumption to vary by
#' treatment status with the first entry being for control and the
#' second for treated.
#' @param c_se A numeric vector of positive values that allow for the
#' bounds to be estimated under a relaxation of `c_se[i]` standard
#' errors of the population constraints. Bounds are returned for the
#' smallest value of `c_se` where bounds are feasible. By default,
#' this is set to a grid of `seq(0, 3, 0.1)`.
#'
#' @return A list object containing bounds.
#'
#' @examples
#' x <- "alfa,bravo,charlie,delta"
#' @export
prepost_bounds <- function(formula, data,  moderator,  prepost,
                           sims = 1000, conf_level = 0.95,
                           moderator_mono = FALSE, outcome_mono = FALSE,
                           stable_mod = FALSE,
                           c_se = NULL,
                           progress = TRUE) {

  outcome = all.vars(formula)[1]
  treat = all.vars(formula)[2]
  moderator = all.vars(moderator)[1]
  order_var = all.vars(prepost)[1]


  Y <- data[, outcome]
  D <- data[, moderator]
  T <- data[, treat]
  Z <- data[, order_var]



  post_est <- mean(Y[T == 1 & D == 1 & Z == 1]) -
    mean(Y[T == 0 & D == 1 & Z == 1]) -
    mean(Y[T == 1 & D == 0 & Z == 1]) +
    mean(Y[T == 0 & D == 0 & Z == 1])
  pre_est <- mean(Y[T == 1 & D == 1 & Z == 0]) -
    mean(Y[T == 0 & D == 1 & Z == 0]) -
    mean(Y[T == 1 & D == 0 & Z == 0]) +
    mean(Y[T == 0 & D == 0 & Z == 0])


  obs <- compute_strata_probs(Y = Y, D = D, T = T, Z = Z)
  mm <- if (moderator_mono) obs$MM else NULL
  om <- if (outcome_mono) obs$OM else NULL

  bounds <- prepost_bounds_core(P = obs$P, V = obs$V, Q = obs$Q, W = obs$W,
                                P_den = obs$P_den, V_den = obs$V_den,
                                outcome_mono = om,
                                moderator_mono = mm,
                                stable_mod = stable_mod)

  out <- list()
  out$Q <- bounds$Q
  out$upper <- bounds$upper
  out$lower <- bounds$lower
  out$c_se <- bounds$c_se

  boot_lo <- boot_hi <- rep(NA, sims)
  if (progress) cat("Bootstrap running...\n")
  empty_count <- 0
  for (b in 1:sims) {
    if ((100 * b / sims) %% 10 == 0 & progress) {
      pb <-  c(rep("==", times = (10 * b / sims)),
               rep("  ", times = 10 - (10 * b / sims)))
      pb <- paste0(c("0% |", pb, "| 100%\n"), collapse = "")

      cat(pb)
    }

    good_draw <- FALSE
    count <- 0
    while (!good_draw) {
      star <- sample(seq_len(nrow(data)), size = nrow(data), replace = TRUE)
      Ys <- Y[star]
      Ds <- D[star]
      Ts <- T[star]
      Zs <- Z[star]
      ostar <- compute_strata_probs(Y = Ys, D = Ds,
                                    T = Ts, Z = Zs)
      good_draw <- !any(is.na(c(ostar$P, ostar$V)))
      count <- count + 1L
      if (count == 100L) {
        stop("100 bootstrap resamples with empty cells.")
      }
    }
    if (count > 1) empty_count <- empty_count + 1

    mm <- if (moderator_mono) ostar$MM else NULL
    om <- if (outcome_mono) ostar$OM else NULL
    bounds <- prepost_bounds_core(
      P = ostar$P, V = ostar$V, Q = ostar$Q, W = ostar$W,
      ostar$P_den, ostar$V_den,
      moderator_mono = mm,
      outcome_mono = mm,
      stable_mod = stable_mod,
      c_se = c_se
    )

    boot_hi[b] <- bounds$upper
    boot_lo[b] <- bounds$lower
  }
  if (empty_count > 0) {
    warning(empty_count, " bootstrap replications required resampling due to empty cells.", call. = FALSE)
  }
  boot_sd_lo <- sd(boot_lo, na.rm = TRUE)
  boot_sd_hi <- sd(boot_hi, na.rm = TRUE)
  boot_ci <- imb.man.ci(out$lower, out$upper, boot_sd_lo, boot_sd_hi,
                         N = nrow(data), alpha = conf_level)

  out$ci_lower <- boot_ci[1]
  out$ci_upper <- boot_ci[2]
  out$pre_est <- pre_est
  out$post_est <- post_est
  return(out)
}



#' Run sensitivity analysis for the randomized moderator placement design
#'
#' @inheritParams post_sens
#' @param prepost A one-sided formula with syntax ~ z, where z is the indicator variable for whether the moderator was measured pre- or post-treatment. 
#' @param t_by Numeric indicating the grid spacing for the
#' \eqn{\theta} parameter that restricts what proportion of units have
#' their outcomes affected by the pre vs post-measurement of the
#' moderator.
#' @param outcome_mono A integer or vector of length 2 indicating
#' if the bounds should assume monotonicity of the effect of the
#' post-test on the outcome with `1` indicating that the post-test
#' effect is positive and `-1` indicating that it is negative. The
#' vector of length 2 allows the monotonicity assumption to vary by
#' treatment status with the first entry being for control and the
#' second for treated.
#' @param c_se A numeric vector of positive values that allow for the
#' bounds to be estimated under a relaxation of `c_se[i]` standard
#' errors of the population constraints. Bounds are returned for the
#' smallest value of `c_se` where bounds are feasible. By default,
#' this is set to a grid of `seq(0, 3, 0.1)`.
#'
#' @return A list object containing sensitivity output.
#'
#' @examples
#' x <- "alfa,bravo,charlie,delta"
#' @export
prepost_sens <- function(formula, data, moderator, prepost,
                         g_by, t_by, sims = 1000,
                         conf_level = 0.95, moderator_mono = NULL,
                         outcome_mono = NULL, c_se = NULL,
                         progress = TRUE) {


  outcome = all.vars(formula)[1]
  treat = all.vars(formula)[2]
  moderator = all.vars(moderator)[1]
  order_var = all.vars(prepost)[1]

  Y <- data[, outcome]
  D <- data[, moderator]
  T <- data[, treat]
  Z <- data[, order_var]

  obs <- compute_strata_probs(Y = Y, D = D, T = T, Z = Z)
  no_con <- prepost_bounds_core(P = obs$P, V = obs$V, Q = obs$Q, W = obs$W,
                                P_den = obs$P_den, V_den = obs$V_den)
  rhos <- seq(0, no_con$max_rho, by = g_by)
  thetas <- seq(0, no_con$max_theta, by = t_by)
  sens_out <- list()
  sens_out$rho <- rhos
  sens_out$theta <- thetas
  sens_out$lower <- matrix(NA, nrow = length(rhos), ncol = length(thetas))
  sens_out$upper <- matrix(NA, nrow = length(rhos), ncol = length(thetas))
  for (r in seq_along(rhos)) {
    for (tt in seq_along(thetas)) {
      out <- prepost_bounds_core(P = obs$P, V = obs$V, Q = obs$Q, W = obs$W,
                                 P_den = obs$P_den, V_den = obs$V_den,
                                 rho = rhos[r], theta = thetas[tt],
                                 moderator_mono = moderator_mono,
                                 outcome_mono = outcome_mono,
                                 c_se = c_se)
      sens_out$lower[r, tt] <- out$lower
      sens_out$upper[r, tt] <- out$upper
    }
  }
  sens_lo <- array(NA, c(sims, length(sens_out$rho), length(sens_out$theta)))
  sens_hi <- array(NA, c(sims, length(sens_out$rho), length(sens_out$theta)))
  for (b in 1:sims) {
    if ((100 * b / sims) %% 10 == 0 & progress) {
      pb <-  c(rep("==", times = (10 * b / sims)),
               rep("  ", times = 10 - (10 * b / sims)))
      pb <- paste0(c("0% |", pb, "| 100%\n"), collapse = "")

      cat(pb)
    }

    star <- sample(seq_len(nrow(data)), size = nrow(data), replace = TRUE)
    Ys <- Y[star]
    Ds <- D[star]
    Ts <- T[star]
    Zs <- Z[star]

    ostar <- compute_strata_probs(Y = Ys, D = Ds, T = Ts, Z = Zs)
    for (r in seq_along(rhos)) {
      for (tt in seq_along(thetas)) {
        out <- prepost_bounds_core(
          P = ostar$P, V = ostar$V, Q = ostar$Q, W = ostar$W,
          ostar$P_den, ostar$V_den, rho = rhos[r], theta = thetas[tt],
          moderator_mono = moderator_mono,
          outcome_mono = outcome_mono,
          c_se = c_se
        )
        sens_lo[b, r, tt] <- out$lower
        sens_hi[b, r, tt] <- out$upper
      }
    }
  }
  sd_lo <- apply(sens_lo, c(2, 3), sd, na.rm = TRUE)
  sd_hi <- apply(sens_hi, c(2, 3), sd, na.rm = TRUE)
  sens_ci_lo <- matrix(NA, nrow = length(sens_out$rho),
                       ncol = length(sens_out$theta))
  sens_ci_hi <- matrix(NA, nrow = length(sens_out$rho),
                       ncol = length(sens_out$theta))
  for (i in seq_along(sens_out$rho)) {
    for (tt in seq_along(sens_out$theta)) {
      if (!is.na(sens_out$lower[i, tt])) {
        imbman <- imb.man.ci(sens_out$lower[i, tt], sens_out$upper[i, tt],
                             sd_lo[i, tt], sd_hi[i, tt], N = nrow(data),
                             alpha = conf_level)
        sens_ci_lo[i, tt] <- imbman[1]
        sens_ci_hi[i, tt] <- imbman[2]
      }
    }
  }
  return(list(rho = sens_out$rho, theta = sens_out$theta,
              lower = sens_out$lower,
              upper = sens_out$upper, ci_lower = sens_ci_lo,
              ci_upper = sens_ci_hi))
}


#' @importFrom lpSolve lp
prepost_bounds_core <- function(P, V, Q, W, P_den, V_den,
                                moderator_mono = NULL,
                                outcome_mono = NULL, stable_mod = FALSE,
                                rho = 1, theta = 1,
                                c_se = NULL) {

  ## P is missing y0 and d0 so we add those
  ## P_Y1_D1_T
  P_levs <- list(Y1 = c(1, 0), D1 = c(1, 0), T = c(1, 0))
  P_vals <- expand.grid(P_levs)
  ## psi_Y1_Y0_D1_T_D0
  psi_levs <- list(Y1 = c(1, 0), Y0 = c(1, 0), D1 = c(1, 0), T = c(1, 0),
                   D0 = c(1, 0))
  psis <- expand.grid(psi_levs)
  num_strata <- nrow(psis)
  Peqmat <- matrix(0, nrow = length(P), ncol = num_strata)
  colnames(Peqmat) <- do.call(paste0, psis)
  for (k in seq_along(P)) {
    pos1 <- which(psis$Y1 == P_vals$Y1[k] & psis$D1 == P_vals$D1[k] &
                    psis$T == P_vals$T[k] & psis$D0 == 1)
    pos0 <- which(psis$Y1 == P_vals$Y1[k] & psis$D1 == P_vals$D1[k] &
                    psis$T == P_vals$T[k] & psis$D0 == 0)
    Peqmat[k, pos1] <- Q
    Peqmat[k, pos0] <- 1 - Q
  }
  ## V_Y0_T_D0
  Veqmat <- matrix(0, nrow = length(V), ncol = num_strata)
  V_levs <- list(Y0 = c(1, 0), T = c(1, 0), D0 = c(1, 0))
  V_vals <- expand.grid(V_levs)
  colnames(Veqmat) <- colnames(Peqmat)
  for (k in seq_along(V)) {
    pos <- which(psis$Y0 == V_vals$Y0[k] & psis$T == V_vals$T[k] &
                   psis$D0 == V_vals$D0[k])
    Veqmat[k, pos] <- 1
  }

  ## sum constraints within TD strata
  td_grid <- expand.grid(psi_levs[c("T", "D0")])
  sum_con <- matrix(0, nrow = nrow(td_grid), ncol = num_strata)
  colnames(sum_con) <- colnames(Peqmat)
  for (k in seq_len(nrow(td_grid))) {
    pos <- which(psis$T == td_grid$T[k] & psis$D0 == td_grid$D0[k])
    sum_con[k, pos] <- 1
  }

  ## rho constraints with TX strata
  t_grid <- expand.grid(T = c(1, 0))
  rho_con <- matrix(0, nrow = nrow(t_grid), ncol = num_strata)
  colnames(rho_con) <- colnames(Peqmat)
  for (k in seq_len(nrow(t_grid))) {
    pos01 <- which(psis$D1 == 0 & psis$T == t_grid$T[k] & psis$D0 == 1)
    pos10 <- which(psis$D1 == 1 & psis$T == t_grid$T[k] & psis$D0 == 0)
    rho_con[k, pos10] <- 1 - Q
    rho_con[k, pos01] <- Q
  }

  ## theta constraints within TD0X strata
  theta_con <- matrix(0, nrow = nrow(td_grid), ncol = num_strata)
  colnames(theta_con) <- colnames(Peqmat)
  for (k in seq_len(nrow(td_grid))) {
    pos <- which(psis$Y1 != psis$Y0 & psis$T == td_grid$T[k] &
                   psis$D0 == td_grid$D0[k])
    theta_con[k, pos] <- 1
  }


  f_obj <- rep(0, times = num_strata)
  f_obj[psis$Y1 == 1 & psis$T == 1 & psis$D0 == 1] <- 1
  f_obj[psis$Y1 == 1 & psis$T == 0 & psis$D0 == 1] <- -1
  f_obj[psis$Y1 == 1 & psis$T == 1 & psis$D0 == 0] <- -1
  f_obj[psis$Y1 == 1 & psis$T == 0 & psis$D0 == 0] <- 1


  zeros <- numeric(0)
  if (stable_mod) {
    sm_pos <- which(psis$D0 != psis$D1 & psis$T == 0)
    zeros <- c(zeros, sm_pos)
    ## sm_con <- matrix(0, nrow = 1, ncol = num_strata)
    ## sm_con[which(psis$T == 0 & psis$D1 == 1 & psis$D0 == 0)] <- 1
    ## sm_con[which(psis$T == 0 & psis$D1 == 0 & psis$D0 == 1)] <- -1
    ## sm_dir <- "=="
  }
  if (!is.null(moderator_mono)) {
    mm_pos <- vector("list", length(moderator_mono))
    for (j in 1:nrow(t_grid)) {
      if (moderator_mono[j] == 1) {
        mm_pos[[j]] <-  which(psis$D0 == 1 & psis$D1 == 0 &
                                psis$T == t_grid$T[j])
      } else {
        mm_pos[[j]] <-  which(psis$D0 == 0 & psis$D1 == 1 &
                                psis$T == t_grid$T[j])
      }
    }
    zeros <- c(zeros, unlist(mm_pos))
  }
  if (!is.null(outcome_mono)) {
    om_pos <- vector("list", length(outcome_mono))
    for (j in 1:nrow(t_grid)) {
      if (outcome_mono[j] == 1) {
        om_pos[[j]] <-  which(psis$Y0 == 1 & psis$Y1 == 0 &
                                psis$T == t_grid$T[j])
      } else {
        om_pos[[j]] <-  which(psis$Y0 == 0 & psis$Y1 == 1 &
                                psis$T == t_grid$T[j])
      }
    }
    zeros <- c(zeros, unlist(om_pos))
  }
  zeros <- unique(zeros)

  pos_con <- diag(nrow = num_strata)
  pos_dirs <- ifelse(seq_len(ncol(pos_con)) %in% zeros, "==", ">=")

  Pse <- sqrt(P * (1 - P) / P_den)
  Vse <- sqrt(V * (1 - V) / V_den)
  if (is.null(c_se)) {
    c_se <- seq(0, 3, by = 0.01)
  }
  for (i in seq_along(c_se)) {
    f_dir <- c(
      rep("<=", nrow(Peqmat) + nrow(Veqmat)),
      rep(">=", nrow(Peqmat) + nrow(Veqmat)),
      rep("==", nrow(sum_con)),
      pos_dirs,
      rep("<=", times = nrow(rho_con) + nrow(theta_con))
    )

    if (sum(V_den) == 0) {
      V_hi <- rep(1, times = length(V))
      V_lo <- rep(0, times = length(V))
    } else {
      V_hi <- V + c_se[i] * Vse
      V_lo <- V - c_se[i] * Vse
    }

    b <- c(
      P + c_se[i] * Pse, V_hi,
      P - c_se[i] * Pse, V_lo,
      rep(1, times = nrow(sum_con)),
      rep(0, times = num_strata),
      rep(rho, times = nrow(rho_con)),
      rep(theta, times = nrow(theta_con))
    )

    A <- rbind(
      Peqmat, Veqmat, Peqmat, Veqmat, sum_con,
      pos_con, rho_con, theta_con
    )

    u_res <- lpSolve::lp("max", f_obj, A, f_dir, b)
    l_res <- lpSolve::lp("min", f_obj, A, f_dir, b)
    lower <- ifelse(l_res$status == 0, l_res$objval, NA)
    upper <- ifelse(u_res$status == 0, u_res$objval, NA)

    ## grab each CATE here
    d1_pos <- which(psis$D0 == 1)
    d0_pos <- which(psis$D0 == 0)
    cate_u_1 <- sum(u_res$solution[d1_pos] * f_obj[d1_pos])
    cate_u_0 <- sum(u_res$solution[d0_pos] * f_obj[d0_pos])
    cate_l_1 <- sum(l_res$solution[d1_pos] * f_obj[d1_pos])
    cate_l_0 <- sum(l_res$solution[d0_pos] * f_obj[d0_pos])

    ## u_res <- Rglpk::Rglpk_solve_LP(f_obj, A, f_dir, b, max = TRUE)
    ## l_res <- Rglpk::Rglpk_solve_LP(f_obj, A, f_dir, b, max = FALSE)
    ## lower <- ifelse(l_res$status == 0, l_res$optimum, NA)
    ## upper <- ifelse(u_res$status == 0, u_res$optimum, NA)
    if (!is.na(lower) & !is.na(upper)) break
  }

  ## calculate the largest value of rho/theta at bounds
  if (!is.na(lower) & !is.na(upper)) {
    max_rho <- max_theta <- 0
    mm_pos100 <- which(psis$D1 == 1 & psis$T == 0 &
                         psis$D0 == 0)
    mm_pos001 <- which(psis$D1 == 0 & psis$T == 0 &
                         psis$D0 == 1)
    mm_pos110 <- which(psis$D1 == 1 & psis$T == 1 &
                         psis$D0 == 0)
    mm_pos011 <- which(psis$D1 == 0 & psis$T == 1 &
                         psis$D0 == 1)
    rho1_up <- sum(u_res$solution[mm_pos110] * (1 - Q) +
                     u_res$solution[mm_pos011] * Q)
    rho0_up <- sum(u_res$solution[mm_pos100] * (1 - Q) +
                     u_res$solution[mm_pos001] * Q)
    rho1_lo <- sum(l_res$solution[mm_pos110] * (1 - Q) +
                     l_res$solution[mm_pos011] * Q)
    rho0_lo <- sum(l_res$solution[mm_pos100] * (1 - Q) +
                     l_res$solution[mm_pos001] * Q)
    max_rho <- max(max_rho, max(rho1_up, rho0_up, rho1_lo, rho0_lo))

    om_pos11 <- which(psis$Y1 != psis$Y0 & psis$T == 1 & psis$D0 == 1)
    om_pos01 <- which(psis$Y1 != psis$Y0 & psis$T == 0 & psis$D0 == 1)
    om_pos10 <- which(psis$Y1 != psis$Y0 & psis$T == 1 & psis$D0 == 0)
    om_pos00 <- which(psis$Y1 != psis$Y0 & psis$T == 0 & psis$D0 == 0)

    thetas <- rep(0, times = 8)
    thetas[1] <- sum(u_res$solution[om_pos11])
    thetas[2] <- sum(u_res$solution[om_pos01])
    thetas[3] <- sum(u_res$solution[om_pos10])
    thetas[4] <- sum(u_res$solution[om_pos00])
    thetas[5] <- sum(l_res$solution[om_pos11])
    thetas[6] <- sum(l_res$solution[om_pos01])
    thetas[7] <- sum(l_res$solution[om_pos10])
    thetas[8] <- sum(l_res$solution[om_pos00])

    max_theta <- max(max_theta, max(thetas))
  } else {
    max_rho <- max_theta <- NA
  }
  return(list(Q = Q, lower = lower, upper = upper, c_se = c_se[i],
              max_rho = max_rho, max_theta = max_theta,
              cate_1 = c(lower = cate_l_1, upper = cate_u_1),
              cate_0 = c(lower = cate_l_0, upper = cate_u_0)))
}




pre_bounds_core <- function(V, outcome_mono = NULL, theta = 1) {

  psi_levs <- list(Y1 = c(1, 0), Y0 = c(1, 0), T = c(1, 0),
                   D0 = c(1, 0))
  psis <- expand.grid(psi_levs)
  num_strata <- nrow(psis)
  ## V_Y0_T_D0_X
  Veqmat <- matrix(0, nrow = length(V), ncol = num_strata)
  V_levs <- list(Y0 = c(1, 0), T = c(1, 0), D0 = c(1, 0))
  V_vals <- expand.grid(V_levs)
  colnames(Veqmat) <- do.call(paste0, psis)
  for (k in seq_along(V)) {
    pos <- which(psis$Y0 == V_vals$Y0[k] & psis$T == V_vals$T[k] &
                   psis$D0 == V_vals$D0[k])
    Veqmat[k, pos] <- 1
  }

  ## sum constraints within TD strata
  td_grid <- expand.grid(psi_levs[c("T", "D0")])
  sum_con <- matrix(0, nrow = nrow(td_grid), ncol = num_strata)
  colnames(sum_con) <- colnames(Veqmat)
  for (k in seq_len(nrow(td_grid))) {
    pos <- which(psis$T == td_grid$T[k] & psis$D0 == td_grid$D0[k])
    sum_con[k, pos] <- 1
  }

  ## theta constraints within TD0X strata
  theta_con <- matrix(0, nrow = nrow(td_grid), ncol = num_strata)
  colnames(theta_con) <- colnames(Veqmat)
  for (k in seq_len(nrow(td_grid))) {
    pos <- which(psis$Y1 != psis$Y0 & psis$T == td_grid$T[k] &
                   psis$D0 == td_grid$D0[k])
    theta_con[k, pos] <- 1
  }

  f_obj <- rep(0, times = num_strata)
  f_obj[which(psis$Y1 == 1 & psis$T == psis$D0)] <- 1
  f_obj[which(psis$Y1 == 1 & psis$T != psis$D0)] <- -1

  zeros <- numeric(0)
  if (!is.null(outcome_mono)) {
    if (length(outcome_mono) == 1) outcome_mono <- rep(outcome_mono, 2)
    if (outcome_mono[1] == 1) {
      om_pos0 <- which(psis$Y1 == 0 & psis$Y0 == 1 & psis$T == 0)
    } else {
      om_pos0 <- which(psis$Y1 == 1 & psis$Y0 == 0 & psis$T == 0)
    }
    if (outcome_mono[2] == 1) {
      om_pos1 <- which(psis$Y1 == 0 & psis$Y0 == 1 & psis$T == 1)
    } else {
      om_pos1 <- which(psis$Y1 == 1 & psis$Y0 == 0 & psis$T == 1)
    }
    zeros <- c(zeros, om_pos0, om_pos1)
  }
  zeros <- unique(zeros)

  pos_con <- diag(nrow = num_strata)
  pos_dirs <- ifelse(seq_len(ncol(pos_con)) %in% zeros, "=", ">=")
  f_dir <- c(rep("=", nrow(Veqmat)),
             rep("=", + nrow(sum_con)),
             pos_dirs,
             rep("<=", times = nrow(theta_con)))
  b <- c(V, rep(1, times = nrow(sum_con)),
         rep(0, times = num_strata),
         rep(theta, times = nrow(theta_con)))
  A <- rbind(Veqmat, sum_con, pos_con, theta_con)
  u_res <- lpSolve::lp("max", f_obj, A, f_dir, b)
  l_res <- lpSolve::lp("min", f_obj, A, f_dir, b)
  lower <- ifelse(l_res$status == 0, l_res$objval, NA)
  upper <- ifelse(u_res$status == 0, u_res$objval, NA)
  return(list(lower = lower, upper = upper))
}


## PRE BOUNDS ----

#' Run pre-treatment bounds.
#'
#' @inheritParams prepost_gibbs
#' @param conf_level A numeric indicating the confidence level for the bootstrap
#'   confidence intervals.
#' @param priming_mono A integer indicating the direction of the priming
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
                        priming_mono = 1L) {


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
  s <- sqrt(p * (1 - p) / q)

  if (priming_mono == 1) {
    s_upper <- s["101"] + s["000"]
    s_lower <- s["100"] + s["001"]
    out$lower <- -(p["100"] + p["001"])
    out$upper <- p["101"] + p["000"]
  } else if (priming_mono == -1) {
    s_upper <- s["001"] + s["100"]
    s_lower <- s["101"] + s["000"]
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
                      priming_mono = 1L) {

  


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
  s <- sqrt(p * (1 - p) / q)  
  
  pre_est <- p["101"] - p["001"] - p["100"] + p["000"]  
  pre_se <- s["101"] + s["001"] + s["100"] + s["000"]
  
  if (priming_mono == 1) {
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
