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
#'   moderator = ~ itaid_bin,
#'  sims = 50
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
  N <- length(Y)

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
  out$lower <- unname(min(bounds$L, na.rm = TRUE))
  out$upper <- unname(max(bounds$U, na.rm = TRUE))
  

  m <- N ## round(N / log(log(N)))
  d_L_vec <- (1 - sqrt(m / N)) * (out$lower - bounds$L)
  d_U_vec <- (1 - sqrt(m / N)) * (out$upper - bounds$U)

  
  if (nondiff) {
    nd_out <- post_bounds_nondiff(p, q)
    out$lower <- nd_out$lower
    out$upper <- nd_out$upper
  }

  boot_lo <- rep(NA, times = sims)
  boot_hi <- rep(NA, times = sims)
  empty_count <- 0
  if (progress) cat("Bootstrap running...\n")
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
      star <- sample(seq_len(N), size = N, replace = TRUE)
      Ys <- Y[star]
      Ds <- D[star]
      Ts <- T[star]
      Zs <- Z[star]

      p <- tapply(Ys, interaction(Ts, Zs, Ds, sep = ""), mean)
      q <- tapply(Ds, interaction(Ts, Zs, sep = ""), mean)
      good_draw <- !any(is.na(p))
      count <- count + 1L
      if (count == 100L) {
        stop("100 bootstrap resamples with empty cells.")
      }
    }
    if (count > 1) empty_count <- empty_count + 1
    
    
    bounds <- post_factory(p, q, moderator_mono, stable_mod)
    ## bounds$L[is.nan(bounds$L)] <- -2
    ## bounds$U[is.nan(bounds$U)] <- 2
    boot_lo[b] <- min(bounds$L, na.rm = TRUE)
    boot_hi[b] <- max(bounds$U, na.rm = TRUE)
    if (nondiff) {
      nd_out <- post_bounds_nondiff(p, q)
      boot_hi[b] <- nd_out$lower
      boot_lo[b] <- nd_out$upper
    }
  }

  boot_sd_lo <- sd(boot_lo, na.rm = TRUE)
  boot_sd_hi <- sd(boot_hi, na.rm = TRUE)
  boot_ci <- imb.man.ci(out$lower, out$upper, boot_sd_lo, boot_sd_hi,
                        N = N, alpha = conf_level)
  out$ci_lower <- max(unname(boot_ci[1]), -2)
  out$ci_upper <- min(unname(boot_ci[2]), 2)

  out$post_est <- post_est
  return(out)
}

post_factory <- function(p, q, moderator_mono, stable_mod) {
  core <- function(Q0) {
    core <- p["111"] * (q["11"] / Q0) - p["011"] * (q["01"] / Q0)
    core <- core - p["110"] * (1 - q["11"]) / (1 - Q0)
    core <- core + p["010"] * (1 - q["01"]) / (1 - Q0)
  }

  ## we remove the NAs because they are from situations when
  ## Q0 is either 0 or 1, which we assume cannot happen so it
  ## must be one of the other points. 
  if (!is.null(moderator_mono)) {
    if (identical(moderator_mono, c(1, 1))) {
      upper <- function(Q0) {
        den <- (Q0 * (1 - Q0))
        core(Q0) + pmin(Q0 - p["111"] * q["11"], 0) / den +
          pmin(p["011"] * q["01"], q["01"] - Q0) / den
      }
      U <- upper(c(
        p["111"] * q["11"],
        q["01"]  * (1 - p["011"]),
        min(q["11"], q["01"])
      ))
      lower <- function(Q0) {
        den <- (Q0 * (1 - Q0))
        core(Q0) - pmin(p["111"] * q["11"], q["11"] - Q0) / den +
          pmax((p["011"] * q["01"] - Q0), 0) / den
      }
      L <- lower(c(
        p["011"] * q["01"], q["11"] * (1 - p["111"]),
        min(q["11"], q["01"])
      ))
    } else if (identical(moderator_mono, c(1, -1))) {
      upper <- function(Q0) {
        den <- Q0 * (1 - Q0)
        core(Q0) - pmax(p["111"] * q["11"] - Q0, 0) / den -
          pmax(p["010"] * (1 - q["01"]) - (1 - Q0), 0) / den
      }
      U <- upper(c(
        p["111"] * q["11"], 1 - p["010"] * (1 - q["01"]),
        q["01"], q["11"]
      ))
      lower <- function(Q0) {
        den <- Q0 * (1 - Q0)
        core(Q0) - pmin(p["111"] * q["11"], q["11"] - Q0) / den -
          pmin(p["010"] * (1 - q["01"]), Q0 - q["01"]) / den
      }
      L <- lower(c(
        q["11"] * (1 - p["111"]),
        p["010"] * (1 - q["01"]) + q["01"],
        q["01"], q["11"]
      ))
    } else if (identical(moderator_mono, c(-1, 1))) {
      upper <- function(Q0){
        den <- Q0 * (1 - Q0)
        core(Q0) + pmin(p["011"] * q["01"], q["01"] - Q0) / den +
           pmin(p["110"] * (1 - q["11"]), Q0 - q["11"]) / den
      }
      U <- upper(c(
        q["01"] * (1 - p["011"]),
        p["110"] * (1 - q["11"]) + q["11"],
        q["01"], q["11"]
      ))
      lower <- function(Q0) {
        den <- Q0 * (1 - Q0)
        core(Q0) + pmax(p["011"] * q["01"] - Q0, 0) / den +
          pmax(p["110"] * (1 - q["11"]) - (1 - Q0), 0) / den

      }
      L <- lower(c(
        p["011"] * q["01"],
        1 - p["110"] * (1 - q["11"]),
        q["01"], q["11"]
      ))
    } else if (identical(moderator_mono, c(-1, -1))) {
      upper <- function(Q0) {
        den <- Q0 * (1 - Q0)
        core(Q0) + pmin(p["110"] * (1 - q["11"]), Q0 - q["11"]) / den -
          pmax(p["010"] * (1 - q["01"]) - (1 - Q0), 0) / den
      }
      U <- upper(c(
        q["11"] + p["110"] * (1 - q["11"]),
        1 - p["010"] * (1 - q["01"]),
        max(q["11"], q["01"])
      ))
      lower <- function(Q0) {
        den <- Q0 * (1 - Q0)
        core(Q0) + pmax(p["110"] * (1 - q["11"]) - (1 - Q0), 0) / den -
          pmin(p["010"] * (1 - q["01"]), Q0 - q["01"]) / den
      }
      L <- lower(c(
        1 - p["110"] * (1 - q["11"]),
        q["01"] + p["010"] * (1 - q["01"]),
        max(q["11"], q["01"])
      ))
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
    U <- upper(c(P1, 1 - P0))
    L <- lower(c(1 - P1, P0))
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
#' @param solver A character indicating what linear programming solver to use:
#'   "Rglpk" (the default) or "lpSolve".
#'
#' @return A list object containing sensitivity output.
#'
#' @examples
#' data(delponte)
#' post_sens(formula = angry_bin ~ t_commonality,
#'   data = delponte,
#'   moderator = ~ itaid_bin,
#'   g_by = 0.1,
#'   sims = 50
#' )
#' @export
post_sens <- function(formula, data,  moderator,
                      g_by, g_max = 1, q_by, sims = 1000, conf_level = 0.95,
                      moderator_mono = NULL, stable_mod = FALSE,
                      progress = TRUE, solver = "Rglpk") {


  # Extract, outcome, moderator, treatment, and covariates variables from formulas

  outcome = all.vars(formula)[1]
  treat = all.vars(formula)[2]
  moderator = all.vars(moderator)[1]


  Y <- data[[outcome]]
  D <- data[[moderator]]
  T <- data[[treat]]
  Z <- rep(1, length(Y))
  N <- length(Y)
  
  obs <- compute_strata_probs(Y = Y, D = D, T = T, Z = Z)
  p <- tapply(Y, interaction(T, Z, D, sep = ""), mean)
  q <- tapply(D, interaction(T, Z, sep = ""), mean)

  
  if (identical(moderator_mono, c(1, 1))) {
    q_max <- min(q["11"], q["01"])
    q_min <- 0
  } else if (identical(moderator_mono, c(-1, -1))) {
    q_min <- max(q["11"], q["01"])
    q_max <- 1
  } else if (identical(moderator_mono, c(-1, 1))) {
    q_min <- q["01"]
    q_max <- q["11"]
  } else if (identical(moderator_mono, c(1, -1))) {
    q_min <- q["11"]
    q_max <- q["01"]
  } else {
    q_min <- 0
    q_max <- 1
  }
  if (stable_mod) q_min <- q_max <- q["01"]

  ## get best bounds 
  Q0 <- seq(q_min, q_max, by = q_by)
  uppers <- rep(NA, length(Q0))
  lowers <- rep(NA, length(Q0))
  g_maxes <- rep(NA, length(Q0))

  for (k in seq_along(Q0)) {
    crit <- post_bounds_core(P = obs$P, Q = Q0[k], 
                             moderator_mono = moderator_mono,
                             stable_mod = stable_mod,
                             estimate_id_region = TRUE,
                             criterion_hat = 0,
                             tau = 0, N = N,
                             fix = NULL)
    bounds <- post_bounds_core(P = obs$P, Q = Q0[k], 
                               moderator_mono = moderator_mono,
                               stable_mod = stable_mod,
                               estimate_id_region = FALSE,
                               criterion_hat = crit$criterion,
                               tau = 0, N = N,
                               fix = NULL)
    uppers[k] <- bounds$upper
    lowers[k] <- bounds$lower
    g_maxes[k] <- bounds$max_rho

  }

  g_cut <- max(g_maxes[which.max(uppers)], g_maxes[which.min(lowers)])
  
  min_gamma <- abs(q["11"] - q["01"])
  gammas <- seq(min_gamma, g_max, by = g_by)
  gam_seq <- which(gammas <= g_cut)
  
  sens_out <- list()
  sens_out$gamma <- gammas
  sens_out$lower <- rep(NA, times = length(gammas))
  sens_out$upper <- rep(NA, times = length(gammas))
  for (g in gam_seq) {
    
    this_q_min <- max(0, q["01"] - gammas[g], q["11"] - gammas[g])
    this_q_min <- this_q_min - 2 * sqrt(this_q_min * (1 - this_q_min) / N)
    this_q_min <- max(0, this_q_min)
    this_q_max <- min(1, gammas[g] + q["01"], gammas[g] + q["11"])
    this_q_max <- this_q_max + 2 * sqrt(this_q_max * (1 - this_q_max) / N)
    this_q_max <- min(1, this_q_max)
    if (missing(q_by)) {
      upper <- optimize(gamma_bounds, interval = c(q_min, q_max), maximum = TRUE,
                        p = p, q = q, gamma = gammas[g], hi = TRUE)
      lower <- optimize(gamma_bounds, interval = c(q_min, q_max), maximum = FALSE,
                        p = p, q = q, gamma = gammas[g], hi = FALSE)
      sens_out$upper[g] <- upper$objective
      sens_out$lower[g] <- lower$objective
    } else {
      Q0 <- seq(max(q_min, this_q_min), min(this_q_max, q_max), by = q_by)
      uppers <- rep(NA, length(Q0))
      lowers <- rep(NA, length(Q0))
      for (k in seq_along(Q0)) {
        crit <- post_bounds_core(obs$P, Q0[k], 
                             moderator_mono = moderator_mono,
                             stable_mod = stable_mod,
                             rho = gammas[g],
                             estimate_id_region = TRUE,
                             criterion_hat = 0,
                             tau = 0.25, N,
                             fix = NULL, solver = solver)
        bounds <- post_bounds_core(obs$P, Q0[k], 
                             moderator_mono = moderator_mono,
                             stable_mod = stable_mod,
                             rho = gammas[g],
                             estimate_id_region = FALSE,
                             criterion_hat = crit$criterion,
                             tau = 0.25, N,
                             fix = NULL, solver = solver)
        uppers[k] <- bounds$upper
        lowers[k] <- bounds$lower
      }
      ## uppers <- gamma_bounds(Q0, p = p, q = q, gamma = gammas[g], hi = TRUE)
      ## lowers <- gamma_bounds(Q0, p = p, q = q, gamma = gammas[g], hi = FALSE)
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
    star <- sample(seq_len(N), size = N, replace = TRUE)
    Ys <- Y[star]
    Ds <- D[star]
    Ts <- T[star]
    Zs <- Z[star]
    ostar <- compute_strata_probs(Y = Ys, D = Ds, T = Ts, Z = Zs)

    ps <- tapply(Ys, interaction(Ts, Zs, Ds, sep = ""), mean)
    qs <- tapply(Ds, interaction(Ts, Zs, sep = ""), mean)

    for (g in gam_seq) {
      this_q_min <- max(0, qs["01"] - gammas[g], qs["11"] - gammas[g])
      this_q_min <- this_q_min - 2 * sqrt(this_q_min * (1 - this_q_min) / N)
      this_q_min <- max(0, this_q_min)
      this_q_max <- min(1, gammas[g] + qs["01"], gammas[g] + qs["11"])
      this_q_max <- this_q_max + 2 * sqrt(this_q_max * (1 - this_q_max) / N)
      this_q_max <- min(1, this_q_max)
      
      if (q_min < q_max) {
        if (missing(q_by)) {
          upper <- optimize(gamma_bounds, interval = c(q_min, q_max), maximum = TRUE,
                            p = ps, q = qs, gamma = gammas[g], hi = TRUE)
          lower <- optimize(gamma_bounds, interval = c(q_min, q_max), maximum = FALSE,
                            p = ps, q = qs, gamma = gammas[g], hi = FALSE)
          sens_hi[b, g] <- upper$objective
          sens_lo[b, g] <- lower$objective
        } else {
          ## if gamma is infeasible, just search over the endpoints
          if (this_q_min > this_q_max) {
            Q0 <- seq(this_q_max, this_q_min, by = q_by)
          } else {
            Q0 <- seq(this_q_min, this_q_max, by = q_by)
          }
          
          uppers <- rep(NA, length(Q0))
          lowers <- rep(NA, length(Q0))
          for (k in seq_along(Q0)) {
            crit <- post_bounds_core(
              ostar$P, Q0[k], 
              moderator_mono = moderator_mono,
              stable_mod = stable_mod,
              rho = gammas[g],
              estimate_id_region = TRUE,
              criterion_hat = 0,
              tau = 0.25, N,
              fix = NULL, solver = solver)
            bounds <- post_bounds_core(
              ostar$P, Q0[k], 
              moderator_mono = moderator_mono,
              stable_mod = stable_mod,
              rho = gammas[g],
              estimate_id_region = FALSE,
              criterion_hat = crit$criterion,
              tau = 0.25, N,
              fix = NULL, solver = solver)
            uppers[k] <- bounds$upper
            lowers[k] <- bounds$lower
          }
          ## uppers <- gamma_bounds(Q0, p = ps, q = qs, gamma = gammas[g], hi = TRUE)
          ## lowers <- gamma_bounds(Q0, p = ps, q = qs, gamma = gammas[g], hi = FALSE)
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

  more_gams <- which(gammas > g_cut)
  sens_out$upper[more_gams] <- sens_out$upper[max(gam_seq)]
  sens_out$lower[more_gams] <- sens_out$lower[max(gam_seq)]
  sens_out$ci_hi[more_gams] <- sens_out$ci_hi[max(gam_seq)]
  sens_out$ci_lo[more_gams] <- sens_out$ci_lo[max(gam_seq)]

  
  out <- data.frame(
    gamma = sens_out$gamma,
    lower = sens_out$lower,
    upper = sens_out$upper,
    ci_lower = sens_out$ci_lo,
    ci_upper = sens_out$ci_hi
  )
  return(out)

}

post_bounds_core <- function(P, Q, 
                             moderator_mono = NULL,
                             stable_mod = FALSE,
                             rho = 1,
                             estimate_id_region = FALSE,
                             criterion_hat = 0,
                             tau = 0.25, N,
                             fix = NULL, solver = "Rglpk") {

  ## P is missing y0 and d0 so we add those
  ## P_Y1_D1_T
  P_levs <- list(Y1 = c(1, 0), D1 = c(1, 0), T = c(1, 0))
  P_vals <- expand.grid(P_levs)
  ## psi_Y1_Y0_D1_T_D0
  psi_levs <- list(Y1 = c(1, 0), Y0 = c(1, 0), D1 = c(1, 0), T = c(1, 0),
                   D0 = c(1, 0))
  psis <- expand.grid(psi_levs)

  ## deterministic constraints: sum to 1 in td_grid
  ## inequality contraints: gamma in t_grid, theta in td_grid  
  td_grid <- expand.grid(psi_levs[c("T", "D0")])
  t_grid <- expand.grid(T = c(1, 0))

  ## we need slack for both theta and gamma constraints
  slack_vars <- nrow(t_grid)
  num_strata <- nrow(psis)
  
  eq_params <- length(P)
  ineq_params <- nrow(t_grid) 
  num_params <- num_strata + 2 * eq_params + 2 * ineq_params + slack_vars    

  ## indices for extra variables
  P_mom <- seq_len(length(P)) + num_strata
  P_abs <- seq_len(length(P)) + max(P_mom)
  rho_mom <- seq_len(nrow(t_grid)) + max(P_abs)
  rho_abs <- seq_len(nrow(t_grid)) + max(rho_mom)
  rho_slack <- seq_len(nrow(t_grid)) + max(rho_abs) 

  ## these will construct the moment conditions
  ## the {P/V}_ind will be the difference between the sum of the psi
  ## values and the corresponding P/V value. 
  Peqmat <- matrix(0, nrow = length(P), ncol = num_params)
  for (k in seq_along(P)) {
    pos1 <- which(psis$Y1 == P_vals$Y1[k] & psis$D1 == P_vals$D1[k] &
                    psis$T == P_vals$T[k] & psis$D0 == 1)
    pos0 <- which(psis$Y1 == P_vals$Y1[k] & psis$D1 == P_vals$D1[k] &
                    psis$T == P_vals$T[k] & psis$D0 == 0)
    Peqmat[k, pos1] <- Q
    Peqmat[k, pos0] <- 1 - Q
    Peqmat[k, P_mom[k]] <- 1
  }

  ## these encode the absolute value transformation.
  ## they will have <= 0 so that the Pa_ind values
  ## will be greater than P_ind or -P_ind (ie, abs(P_ind))
  Pabsmat <- matrix(0, nrow = 2 * length(P), ncol = num_params)
  for (k in seq_along(P)) {
    Pabsmat[k, P_mom[k]] <- -1
    Pabsmat[k, P_abs[k]] <- 1
    Pabsmat[length(P) + k, P_mom[k]] <- 1
    Pabsmat[length(P) + k, P_abs[k]] <- 1

  }
  
  ## sum constraints within TD strata
  sum_con <- matrix(0, nrow = nrow(td_grid), ncol = num_params)
  colnames(sum_con) <- colnames(Peqmat)
  for (k in seq_len(nrow(td_grid))) {
    pos <- which(psis$T == td_grid$T[k] & psis$D0 == td_grid$D0[k])
    sum_con[k, pos] <- 1
  }

  ## rho constraints with TX strata
  rho_con <- matrix(0, nrow = nrow(t_grid), ncol = num_params)
  for (k in seq_len(nrow(t_grid))) {
    pos01 <- which(psis$D1 == 0 & psis$T == t_grid$T[k] & psis$D0 == 1)
    pos10 <- which(psis$D1 == 1 & psis$T == t_grid$T[k] & psis$D0 == 0)
    rho_con[k, pos10] <- 1 - Q
    rho_con[k, pos01] <- Q

    ## slack variables
    rho_con[k, rho_slack[k]] <- 1
    rho_con[k, rho_mom[k]] <- 1
  }

  rho_abs_con <- matrix(0, nrow = 2 * nrow(t_grid), ncol = num_params)
  for (k in seq_len(nrow(t_grid))) {
    rho_abs_con[k, rho_mom[k]] <- -1
    rho_abs_con[k, rho_abs[k]] <- 1
    rho_abs_con[nrow(t_grid) + k, rho_mom[k]] <- 1
    rho_abs_con[nrow(t_grid) + k, rho_abs[k]] <- 1    
  }
  

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
  zeros <- unique(zeros)

  zeros <- unique(zeros)
  zero_con <- matrix(0, nrow = length(zeros), ncol = num_params)
  for (k in seq_along(zeros)) {
    zero_con[k, zeros[k]] <- 1
  }
  zero_dir <- rep("==", times = length(zeros))
  zero_b <- rep(0, times = length(zeros))

  f_dir <- c(
    rep("==", nrow(Peqmat)),
    rep(">=", nrow(Pabsmat)),
    rep("==", nrow(sum_con)),
    zero_dir,
    rep("==", times = nrow(rho_con)),
    rep(">=", times = nrow(rho_abs_con))
  )

  b <- c(
    P,
    rep(0, times = nrow(Pabsmat)),
    rep(1, times = nrow(sum_con)),
    zero_b,
    rep(rho, times = nrow(rho_con)),
    rep(0, times = nrow(rho_abs_con))
    )

  A <- rbind(
    Peqmat, Pabsmat, sum_con,
    zero_con, rho_con, rho_abs_con
  )

  unbounded <- c(P_mom, rho_mom)
  lp_bounds <- list(
    lower = list(ind = unbounded,
                 val = rep(-Inf, length(unbounded))),
    upper = list(ind = unbounded,
                 val = rep(Inf, length(unbounded)))
  )
  
  lower <- NA
  upper <- NA
  cate_u_1 <- NA
  cate_u_0 <- NA
  cate_l_1 <- NA
  cate_l_0 <- NA
  criterion <- NA
  max_rho <- NA
  max_theta <- NA
  solution <- rep(NA, times = num_params)
  solution_L <- rep(NA, times = num_params)
  solution_U <- rep(NA, times = num_params)
  f_obj <- rep(0, times = num_params)
  if (estimate_id_region) {

    if (!is.null(fix)) {
      fix_con <- matrix(0, nrow = 1, ncol = num_params)
      fix_con[which(psis$Y1 == 1 & psis$T == 1 & psis$D0 == 1)] <- 1
      fix_con[which(psis$Y1 == 1 & psis$T == 0 & psis$D0 == 1)] <- -1
      fix_con[which(psis$Y1 == 1 & psis$T == 1 & psis$D0 == 0)] <- -1
      fix_con[which(psis$Y1 == 1 & psis$T == 0 & psis$D0 == 0)] <- 1

      A <- rbind(A, fix_con)
      f_dir <- c(f_dir, "==")
      b <- c(b, fix)
    }
    
    f_obj[P_abs] <- sqrt(N)
    f_obj[rho_abs] <- sqrt(N)

    
    min_crit <- Rglpk::Rglpk_solve_LP(
      f_obj, A, f_dir, b, max = FALSE,
      bounds = lp_bounds,
      control = list(presolve = TRUE)
    )
    criterion <- ifelse(min_crit$status == 0, min_crit$optimum, NA)
    solution <- min_crit$solution
  } else {
    
    f_obj[which(psis$Y1 == 1 & psis$T == 1 & psis$D0 == 1)] <- 1
    f_obj[which(psis$Y1 == 1 & psis$T == 0 & psis$D0 == 1)] <- -1
    f_obj[which(psis$Y1 == 1 & psis$T == 1 & psis$D0 == 0)] <- -1
    f_obj[which(psis$Y1 == 1 & psis$T == 0 & psis$D0 == 0)] <- 1

    crit_con <- matrix(0, nrow = 1, ncol = num_params)
    crit_con[P_abs] <- sqrt(N)
    crit_con[rho_abs] <- sqrt(N)
    A <- rbind(A, crit_con)
    b <- c(b, criterion_hat * (1 + tau))
    f_dir <- c(f_dir, "<=")
    if (all(is.finite(A)) & all(is.finite(b))) {
      ## ## grab each CATE here
      d1_pos <- which(psis$D0 == 1)
      d0_pos <- which(psis$D0 == 0)
            if (solver == "lpSolve") {
        u_res <- lpSolve::lp("max", f_obj, A, f_dir, b)
        l_res <- lpSolve::lp("min", f_obj, A, f_dir, b)
        lower <- ifelse(l_res$status == 0, l_res$objval, NA)
        upper <- ifelse(u_res$status == 0, u_res$objval, NA)

      } else if (solver == "Rglpk") {
        u_res <- Rglpk::Rglpk_solve_LP(f_obj, A, f_dir, b, max = TRUE,
                                       bounds = lp_bounds,
                                       control = list(presolve = TRUE))
        l_res <- Rglpk::Rglpk_solve_LP(f_obj, A, f_dir, b, max = FALSE,
                                       bounds = lp_bounds,
                                       control = list(presolve = TRUE))
        lower <- ifelse(l_res$status == 0, l_res$optimum, NA)
        upper <- ifelse(u_res$status == 0, u_res$optimum, NA)

      }

      cate_u_1 <- sum(u_res$solution[d1_pos] * f_obj[d1_pos])
      cate_u_0 <- sum(u_res$solution[d0_pos] * f_obj[d0_pos])
      cate_l_1 <- sum(l_res$solution[d1_pos] * f_obj[d1_pos])
      cate_l_0 <- sum(l_res$solution[d0_pos] * f_obj[d0_pos])
      criterion <- sum(crit_con * u_res$solution)
      solution_L <- l_res$solution
      solution_U <- u_res$solution
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

    } else {
      max_rho <- max_theta <- NA
    }

  }

  
  return(list(Q = Q, lower = lower, upper = upper,
              max_rho = max_rho, 
              cate_1 = c(lower = cate_l_1, upper = cate_u_1),
              cate_0 = c(lower = cate_l_0, upper = cate_u_0),
              criterion = criterion, solution = solution,
              solution_L = solution_L, solution_U = solution_U))
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
