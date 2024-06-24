



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
                        m = NULL,
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

  m_sample <- sample(seq_len(N), size = m)
  Ym <- Y[m_sample]
  Dm <- D[m_sample]
  Tm <- T[m_sample]
  Zm <- Z[m_sample]
  p_m <- tapply(Ym, interaction(Tm, Zm, Dm, sep = ""), mean)
  q_m <- tapply(Dm, interaction(Tm, Zm, sep = ""), mean)
  bounds_m <- post_factory(p_m, q_m, moderator_mono, stable_mod)
  lower_m <- unname(min(bounds_m$L, na.rm = TRUE))
  upper_m <- unname(max(bounds_m$U, na.rm = TRUE))
  
  if (nondiff) {
    nd_out <- post_bounds_nondiff(p, q)
    out$lower <- nd_out$lower
    out$upper <- nd_out$upper
  }

  boot_lo <- rep(NA, times = sims)
  boot_hi <- rep(NA, times = sims)
  boot_lo_mod <- rep(NA, times = sims)
  boot_hi_mod <- rep(NA, times = sims)
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
    boot_lo_mod[b] <- min(bounds$L + d_L_vec, na.rm = TRUE)
    boot_hi_mod[b] <- max(bounds$U + d_U_vec, na.rm = TRUE)
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
                             fix = NULL)
        bounds <- post_bounds_core(obs$P, Q0[k], 
                             moderator_mono = moderator_mono,
                             stable_mod = stable_mod,
                             rho = gammas[g],
                             estimate_id_region = FALSE,
                             criterion_hat = crit$criterion,
                             tau = 0.25, N,
                             fix = NULL)
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
              fix = NULL)
            bounds <- post_bounds_core(
              ostar$P, Q0[k], 
              moderator_mono = moderator_mono,
              stable_mod = stable_mod,
              rho = gammas[g],
              estimate_id_region = FALSE,
              criterion_hat = crit$criterion,
              tau = 0.25, N,
              fix = NULL)
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
                             fix = NULL) {

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
      u_res <- Rglpk::Rglpk_solve_LP(f_obj, A, f_dir, b, max = TRUE,
                                     bounds = lp_bounds,
                                     control = list(presolve = TRUE))
      l_res <- Rglpk::Rglpk_solve_LP(f_obj, A, f_dir, b, max = FALSE,
                                     bounds = lp_bounds,
                                     control = list(presolve = TRUE))
      lower <- ifelse(l_res$status == 0, l_res$optimum, NA)
      upper <- ifelse(u_res$status == 0, u_res$optimum, NA)
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



post_bounds_inference <- function(P, Q, P_star, Q_star,
                                moderator_mono = NULL,
                                stable_mod = FALSE,
                                rho = 1,
                                estimate_id_region = FALSE,
                                criterion_hat = 0,
                                tau = 0.25, N,
                                fix = NULL) {

  ## P is missing y0 and d0 so we add those
  ## P_Y1_D1_T
  P_levs <- list(Y1 = c(1, 0), D1 = c(1, 0), T = c(1, 0))
  P_vals <- expand.grid(P_levs)
  ## psi_Y1_Y0_D1_T_D0
  psi_levs <- list(Y1 = c(1, 0), Y0 = c(1, 0), D1 = c(1, 0), T = c(1, 0),
                   D0 = c(1, 0))
  psis <- expand.grid(psi_levs)

  P_diff <- P_star - P
  Q_diff <- Q_star - Q
  V_diff <- V_star - V

  ## deterministic constraints: sum to 1 in td_grid
  ## inequality contraints: gamma in t_grid, theta in td_grid  
  td_grid <- expand.grid(psi_levs[c("T", "D0")])
  t_grid <- expand.grid(T = c(1, 0))

  ## we need slack for both theta and gamma constraints
  slack_vars <- nrow(t_grid)
  num_strata <- nrow(psis)
  
  eq_params <- length(P)
  ineq_params <- nrow(t_grid) 
  num_params <- 2 * num_strata + 5 * eq_params + 5 * ineq_params + 2 * slack_vars    

  ## indices for extra variables
  H_par <- seq_len(num_strata) + num_strata
  P_mom <- seq_len(length(P)) + max(H_par)
  P_abs <- seq_len(length(P)) + max(P_mom)
  rho_mom <- seq_len(nrow(t_grid)) + max(P_abs)
  rho_abs <- seq_len(nrow(t_grid)) + max(rho_mom)
  rho_slack <- seq_len(nrow(t_grid)) + max(rho_abs) 
  P_sam <- seq_len(length(P)) + max(rho_slack) 
  P_pos <- seq_len(length(P)) + max(P_sam)
  P_neg <- seq_len(length(P)) + max(P_pos)
  rho_sam <- seq_len(nrow(t_grid)) + max(P_neg) 
  rho_pos <- seq_len(nrow(t_grid)) + max(rho_sam)
  rho_neg <- seq_len(nrow(t_grid)) + max(rho_pos)
  rho_h_slack <- seq_len(nrow(t_grid)) + max(rho_neg)
  
  ## these will construct the moment conditions
  ## the {P/V}_ind will be the difference between the sum of the psi
  ## values and the corresponding P/V value. 
  Peqmat <- matrix(0, nrow = length(P), ncol = num_params)
  for (k in seq_along(P)) {
    pos1 <- which(psis$Y1 == P_vals$Y1[k] & psis$D1 == P_vals$D1[k] &
                    psis$T == P_vals$T[k] & psis$D0 == 1)
    pos0 <- which(psis$Y1 == P_vals$Y1[k] & psis$D1 == P_vals$D1[k] &
                    psis$T == P_vals$T[k] & psis$D0 == 0)
    Peqmat[k, pos1] <- sqrt(N) * Q_diff
    Peqmat[k, pos0] <- -sqrt(N) * Q_diff
    Peqmat[k, P_mom[k]] <- 1

    ## H variables
    Peqmat[k, pos1 + num_strata] <- Q
    Peqmat[k, pos0 + num_strata] <- 1 - Q
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
  ## these encode the absolute value transformation.
  ## they will have <= 0 so that the Pa_ind values
  ## will be greater than P_ind or -P_ind (ie, abs(P_ind))
  Psammat <- matrix(0, nrow = 2 * length(P), ncol = num_params)
  for (k in seq_along(P)) {
    pos1 <- which(psis$Y1 == P_vals$Y1[k] & psis$D1 == P_vals$D1[k] &
                    psis$T == P_vals$T[k] & psis$D0 == 1)
    pos0 <- which(psis$Y1 == P_vals$Y1[k] & psis$D1 == P_vals$D1[k] &
                    psis$T == P_vals$T[k] & psis$D0 == 0)

    Psammat[k, pos1] <- Q
    Psammat[k, pos0] <- 1 - Q
    Psammat[k, P_sam[k]] <- 1
    Psammat[length(P) + k, P_sam[k]] <- 1
    Psammat[length(P) + k, P_pos[k]] <- -1
    Psammat[length(P) + k, P_neg[k]] <- 1
  }
  P_sam_dir <- rep("==", times = nrow(Psammat))
  P_sam_b <- c(P, rep(0, times = length(P)))
  
  
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
    rho_con[k, pos10] <- -sqrt(N) * Q_diff
    rho_con[k, pos01] <- sqrt(N) * Q_diff
    rho_con[k, pos10 + num_strata] <- 1 - Q
    rho_con[k, pos01 + num_strata] <- Q

    ## slack variables
    rho_con[k, rho_h_slack[k]] <- 1
    rho_con[k, rho_mom[k]] <- 1
  }

  rho_abs_con <- matrix(0, nrow = 2 * nrow(t_grid), ncol = num_params)
  for (k in seq_len(nrow(t_grid))) {
    rho_abs_con[k, rho_mom[k]] <- -1
    rho_abs_con[k, rho_abs[k]] <- 1
    rho_abs_con[nrow(t_grid) + k, rho_mom[k]] <- 1
    rho_abs_con[nrow(t_grid) + k, rho_abs[k]] <- 1    
  }
  
  rho_sam_con <- matrix(0, nrow = 2 * nrow(t_grid), ncol = num_params)
  for (k in seq_len(nrow(t_grid))) {
    pos01 <- which(psis$D1 == 0 & psis$T == t_grid$T[k] & psis$D0 == 1)
    pos10 <- which(psis$D1 == 1 & psis$T == t_grid$T[k] & psis$D0 == 0)
    rho_sam_con[k, pos10] <- 1 - Q
    rho_sam_con[k, pos01] <- Q
    rho_sam_con[k, rho_sam[k]] <- 1
    rho_sam_con[k, rho_slack[k]] <- 1
    rho_sam_con[nrow(t_grid) + k, rho_sam[k]] <- 1
    rho_sam_con[nrow(t_grid) + k, rho_pos[k]] <- -1
    rho_sam_con[nrow(t_grid) + k, rho_neg[k]] <- 1
  }

  rho_sam_dir <- rep("==", times = 2 * nrow(t_grid))
  rho_sam_b <- c(rep(rho, times = nrow(t_grid)), rep(0, times = nrow(t_grid)))



  crit_con <- matrix(0, nrow = 1, ncol = num_params)
  crit_con[c(P_pos, P_neg)] <- sqrt(N)
  crit_con[c(rho_pos, rho_neg)] <- sqrt(N)
  crit_dir <- "<="
  crit_b <- criterion_hat * (1 + tau)
  
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
  zero_con <- matrix(0, nrow = 2 * length(zeros), ncol = num_params)
  for (k in seq_along(zeros)) {
    zero_con[k, zeros[k]] <- 1
    zero_con[length(zeros) + k, H_par[zeros[k]]] <- 1
  }
  zero_dir <- rep("==", times = 2 * length(zeros))
  zero_b <- rep(0, times = 2 * length(zeros))

  h_cons <- matrix(0, nrow = 2 * num_strata, ncol = num_params)
  for (k in seq_len(num_strata)) {
    h_cons[k, k] <- 1
    h_cons[k, H_par[k]] <- 1 / sqrt(N)
    h_cons[num_strata + k, k] <- 1
    h_cons[num_strata + k, H_par[k]] <- 1 / sqrt(N)
  }
  rho_slack_con <- matrix(0, nrow = nrow(t_grid), ncol = num_params)
  for (k in seq_len(nrow(t_grid))) {
    rho_slack_con[k, rho_slack[k]] <- 1
    rho_slack_con[k, rho_h_slack[k]] <- 1 / sqrt(N)
  }
  rho_slack_b <- rep(0, times = nrow(t_grid))
  rho_slack_dirs <- rep(">=", times = nrow(t_grid))
  

  h_sum <- matrix(0, nrow = 1, ncol = num_params)
  h_sum[H_par] <- 1

  f_dir <- c(
    rep("==", nrow(Peqmat)),
    rep(">=", nrow(Pabsmat)),
    rep("==", nrow(sum_con)),
    zero_dir,
    rep("==", times = nrow(rho_con)),
    rep(">=", times = nrow(rho_abs_con)),
    rep("<=", times = num_strata),
    rep(">=", times = num_strata),
    P_sam_dir,
    rho_sam_dir,
    crit_dir,
    rho_slack_dirs,
    "=="
  )

  b <- c(
    sqrt(N) * P_diff,
    rep(0, times = nrow(Pabsmat)),
    rep(1, times = nrow(sum_con)),
    zero_b,
    rep(rho, times = nrow(rho_con)),
    rep(0, times = nrow(rho_abs_con)),
    rep(1, times = num_strata),
    rep(0, times = num_strata),
    P_sam_b,
    rho_sam_b,
    crit_b,
    rho_slack_b,
    0
    )

  A <- rbind(
    Peqmat, Pabsmat,
    sum_con,
    zero_con,
    rho_con, rho_abs_con,
    h_cons,
    Psammat,
    rho_sam_con, 
    crit_con,
    rho_slack_con,
    h_sum
  )

  unbounded <- c(H_par, P_mom, rho_mom,
                 P_sam, rho_sam, rho_h_slack)
  lp_bounds <- list(
    lower = list(ind = unbounded,
                 val = rep(-Inf, length(unbounded))),
    upper = list(ind = unbounded,
                 val = rep(Inf, length(unbounded)))
  )
  
  criterion <- NA
  if (!is.null(fix)) {
    fix_con <- matrix(0, nrow = 1, ncol = num_params)
    fix_con[which(psis$Y1 == 1 & psis$T == 1 & psis$D0 == 1)] <- 1
    fix_con[which(psis$Y1 == 1 & psis$T == 0 & psis$D0 == 1)] <- -1
    fix_con[which(psis$Y1 == 1 & psis$T == 1 & psis$D0 == 0)] <- -1
    fix_con[which(psis$Y1 == 1 & psis$T == 0 & psis$D0 == 0)] <- 1

    fix_h_con <- matrix(0, nrow = 1, ncol = num_params)
    fix_h_con[which(psis$Y1 == 1 & psis$T == 1 & psis$D0 == 1) + num_strata] <- 1
    fix_h_con[which(psis$Y1 == 1 & psis$T == 0 & psis$D0 == 1) + num_strata] <- -1
    fix_h_con[which(psis$Y1 == 1 & psis$T == 1 & psis$D0 == 0) + num_strata] <- -1
    fix_h_con[which(psis$Y1 == 1 & psis$T == 0 & psis$D0 == 0) + num_strata] <- 1

    A <- rbind(A, fix_con, fix_con, fix_h_con)
    f_dir <- c(f_dir, "<=", ">=", "==")
    b <- c(b, fix +0, fix - 0, 0)
  }

  f_obj <- rep(0, times = num_params)
  ## f_obj[c(P_pos, P_neg)] <- sqrt(N)
  ## f_obj[c(rho_pos, rho_neg)] <- sqrt(N)

  f_obj[P_abs] <- 1
  f_obj[rho_abs] <- 1
  
  ## browser()
  min_crit <- Rglpk::Rglpk_solve_LP(
    f_obj, A, f_dir, b, max = FALSE,
    bounds = lp_bounds,
    control = list(presolve = TRUE, tm_limit = 10000)
    )
  ## if (is.na(criterion))  browser()
  criterion <- ifelse(min_crit$status == 0, min_crit$optimum, NA)

  
  return(criterion)
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
#' @param prepost A one-sided formula with syntax ~ z, where z is the indicator
#'   variable for whether the moderator was measured pre- or post-treatment.
#' @param sims An integer specifying the number of simulations for the
#'   sensitivity analysis.
#' @param conf_level A numeric specifying level for the confidence intervals.
#' @param outcome_mono A integer or vector of length 2 indicating if the bounds
#'   should assume monotonicity of the effect of the post-test on the outcome
#'   with `1` indicating that the post-test effect is positive and `-1`
#'   indicating that it is negative. The vector of length 2 allows the
#'   monotonicity assumption to vary by treatment status with the first entry
#'   being for control and the second for treated.
#' @param tau A numeric indicating how close the the moment conditions of the
#'   estimated bounds have to be from the minimum values in the sample. This
#'   allows us to obtain bounds and confidence intevals even when the
#'   assumptions are slightly violated due to sampling.
#'
#' @return A list object containing bounds.
#'
#' @examples
#' x <- "alfa,bravo,charlie,delta"
#' @export
prepost_bounds <- function(formula, data,  moderator,  prepost,
                           sims = 1000, conf_level = 0.95,
                           moderator_mono = NULL, outcome_mono = NULL,
                           stable_mod = FALSE,
                           tau = 0.25,
                           progress = TRUE) {

  outcome = all.vars(formula)[1]
  treat = all.vars(formula)[2]
  moderator = all.vars(moderator)[1]
  order_var = all.vars(prepost)[1]


  Y <- data[, outcome]
  D <- data[, moderator]
  T <- data[, treat]
  Z <- data[, order_var]
  N <- length(Y)


  post_est <- mean(Y[T == 1 & D == 1 & Z == 1]) -
    mean(Y[T == 0 & D == 1 & Z == 1]) -
    mean(Y[T == 1 & D == 0 & Z == 1]) +
    mean(Y[T == 0 & D == 0 & Z == 1])
  pre_est <- mean(Y[T == 1 & D == 1 & Z == 0]) -
    mean(Y[T == 0 & D == 1 & Z == 0]) -
    mean(Y[T == 1 & D == 0 & Z == 0]) +
    mean(Y[T == 0 & D == 0 & Z == 0])


  obs <- compute_strata_probs(Y = Y, D = D, T = T, Z = Z)
  ##mm <- if (moderator_mono) obs$MM else NULL
  ##om <- if (outcome_mono) obs$OM else NULL
  om <- outcome_mono
  mm <- moderator_mono

  crit <- prepost_bounds_core(P = obs$P, V = obs$V, Q = obs$Q, W = obs$W,
                                P_den = obs$P_den, V_den = obs$V_den,
                                outcome_mono = om,
                               moderator_mono = mm,
                                stable_mod = stable_mod,
                               estimate_id_region = TRUE, N = N)
  bounds <- prepost_bounds_core(P = obs$P, V = obs$V, Q = obs$Q, W = obs$W,
                                P_den = obs$P_den, V_den = obs$V_den,
                                outcome_mono = om,
                               moderator_mono = mm,
                                stable_mod = stable_mod,
                               criterion_hat = crit$criterion, tau = tau, N = N)

  out <- list()
  out$Q <- bounds$Q
  out$upper <- bounds$upper
  out$lower <- bounds$lower
  out$solution_L <- bounds$solution_L
  out$solution_U <- bounds$solution_U
  out$deviation <- crit$criterion

  ## code to do CNS confidence intervals. Low coverage, but might revisit.
  cns <- FALSE
  if (cns) {
    if (progress) cat("Finding upper CI via CNS...\n")
    upper_ci <- find_endpoint(out, obs, Y, D, T, Z, N, rho = 1, theta = 1, om,
                              mm, stable_mod, sims, upper = TRUE,
                              conf_level, tau, progress)  
    if (progress) cat("Finding lower CI via CNS...\n")
    lower_ci <- find_endpoint(out, obs, Y, D, T, Z, N, rho = 1, theta = 1, om,
                              mm, stable_mod, sims, upper = FALSE,
                              conf_level, tau, progress)
    if (max(lower_ci$p_values) > 1 - conf_level) {    
      out$ci_lower_cns <- min(setdiff(lower_ci$point_list, lower_ci$reject_list))
    } else {
      out$ci_lower_cns <- out$lower
    }
    if (max(upper_ci$p_values) > 1 - conf_level) {
      out$ci_upper_cns <- max(setdiff(upper_ci$point_list, upper_ci$reject_list))
    } else {
      out$ci_upper_cns <- out$upper
    }
  }
  
  boot_lo <- boot_hi <- ms_test <- rep(NA, sims)
  
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
      star <- sample(seq_len(N), size = N, replace = TRUE)
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

    ## mm <- if (moderator_mono) ostar$MM else NULL
    ## om <- if (outcome_mono) ostar$OM else NULL
    crit_star <- prepost_bounds_core(
      P = ostar$P, V = ostar$V, Q = ostar$Q, W = ostar$W,
      P_den = ostar$P_den, V_den = ostar$V_den,
      outcome_mono = om,
      moderator_mono = mm,
      stable_mod = stable_mod,
      estimate_id_region = TRUE, N = N)
    bounds <- prepost_bounds_core(
      P = ostar$P, V = ostar$V, Q = ostar$Q, W = ostar$W,
      P_den = ostar$P_den, V_den = ostar$V_den,
      outcome_mono = om,
      moderator_mono = mm,
      stable_mod = stable_mod,
      criterion_hat = crit_star$criterion, tau = tau, N = N)
    q_ms <- prepost_bounds_inference(
      P = obs$P, V = obs$V, Q = obs$Q,
      P_star = ostar$P, V_star = ostar$V, Q_star = ostar$Q,
      outcome_mono = om,
      moderator_mono = mm,
      stable_mod = stable_mod,
      criterion_hat = crit$criterion, N = N,
      tau= tau)
    ms_test[b] <- q_ms

    boot_hi[b] <- bounds$upper
    boot_lo[b] <- bounds$lower
  }
  if (empty_count > 0) {
    warning(empty_count, " bootstrap replications required resampling due to empty cells.", call. = FALSE)
  }
  out$ms_p <- mean(crit$criterion <= ms_test + 1e-6)
  boot_sd_lo <- sd(boot_lo, na.rm = TRUE)
  boot_sd_hi <- sd(boot_hi, na.rm = TRUE)
  boot_ci <- imb.man.ci(out$lower, out$upper, boot_sd_lo, boot_sd_hi,
                         N = N, alpha = conf_level)
  out$ci_lower <- boot_ci[1]
  out$ci_upper <- boot_ci[2]

  out$pre_est <- pre_est
  out$post_est <- post_est
  return(out)
}



find_endpoint <- function(out, obs, Y, D, T, Z, N, rho, theta,
                          outcome_mono,
                          moderator_mono, stable_mod, sims, upper = TRUE,
                          conf_level = 0.95, tau, progress = TRUE) {

  ## upper CI
  if (upper) {
    hi <- 2
    lo <- out$upper
  } else {
    hi <- out$lower
    lo <- -2
  }
  
  reject_list <- c()
  point_list <- c()
  p_values <- c()

  while ((hi - lo) > 0.01) {
    test <- (hi + lo) / 2
    if (progress) cat("Test point: ", test, "\t")
    q_t <- prepost_bounds_core(
      P = obs$P, V = obs$V, Q = obs$Q, W = obs$W,
      P_den = obs$P_den, V_den = obs$V_den,
      rho = rho, theta = theta,
      outcome_mono = outcome_mono,
      moderator_mono = moderator_mono,
      stable_mod = stable_mod,
      estimate_id_region = TRUE, N = N,
      fix = test
      )
    if (progress)  cat("Obs min: ", q_t$criterion, "\t")
    if (q_t$criterion < 1e-7) {
      this_reject <- FALSE
      p_values <- c(p_values, 1)
    } else {
      bad_mm_count <- rep(0, sims)
      bad_om_count <- rep(0, sims)

      boot_crit <- rep(NA, sims)
      for (b in 1:sims) {
        ## if ((100 * b / sims) %% 10 == 0) {
        ##   pb <-  c(rep("==", times = (10 * b / sims)),
        ##            rep("  ", times = 10 - (10 * b / sims)))
        ##   pb <- paste0(c("0% |", pb, "| 100%\n"), collapse = "")
          
        ##   cat(pb)
        ## } 
        ## cat(b, "\n")
        good_draw <- FALSE
        count <- 0
        empty_count <- 0
        while (!good_draw) {
          star <- sample(seq_len(N), size = N, replace = TRUE)
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
        q_bs <- NA
        tau_use <- tau
        if (!identical(obs$MM, ostar$MM)) bad_mm_count[b] <- 1
        if (!identical(obs$OM, ostar$OM)) bad_om_count[b] <- 1
        while (is.na(q_bs) & tau_use < 1) {
          
          q_bs <- prepost_bounds_inference(
            P = obs$P, V = obs$V, Q = obs$Q,
            P_star = ostar$P, V_star = ostar$V, Q_star = ostar$Q,
            rho = rho, theta = theta,
            outcome_mono = outcome_mono,
            moderator_mono = moderator_mono,
            stable_mod = stable_mod,
            fix = test,
            tau = tau,
            criterion_hat = q_t$criterion, N = N)
          
          if (is.na(q_bs)) {
            tau_use <- tau_use + 0.05
          }

        }
        boot_crit[b] <- q_bs
        
      }
      
      p_values <- c(p_values, mean(q_t$criterion <= boot_crit + 1e-6, na.rm = TRUE))
      this_reject <- mean(q_t$criterion <= boot_crit + 1e-6, na.rm = TRUE) <= (1- conf_level) + 1e-6
      
    } 
    point_list <- c(point_list, test)

    ## if (any(is.nan(p_values))) browser()
    if (progress)  cat("p-value: ", p_values[length(p_values)], "\n")
    
    if (this_reject) {
      if (upper) {
        hi <- test
      } else {
        lo <- test
      }
      reject_list <- c(reject_list, test)
    } else {
      if (upper) {
        lo <- test
      } else {
        hi <- test
      }
      
    }

  }
  return(list(point_list = point_list, p_values = p_values, reject_list = reject_list))


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
#' @param tau A numeric indicating how close the the moment conditions of the
#'   estimated bounds have to be from the minimum values in the sample. This
#'   allows us to obtain bounds and confidence intevals even when the
#'   assumptions are slightly violated due to sampling.
#'
#' @return A list object containing sensitivity output.
#'
#' @examples
#' x <- "alfa,bravo,charlie,delta"
#' @export
prepost_sens <- function(formula, data, moderator, prepost,
                         g_by, g_at, t_by, t_at, sims = 1000, stable_mod = FALSE,
                         conf_level = 0.95, moderator_mono = NULL,
                         outcome_mono = NULL, tau = 0.25,
                         progress = TRUE) {


  outcome = all.vars(formula)[1]
  treat = all.vars(formula)[2]
  moderator = all.vars(moderator)[1]
  order_var = all.vars(prepost)[1]

  Y <- data[, outcome]
  D <- data[, moderator]
  T <- data[, treat]
  Z <- data[, order_var]
  N <- length(Y)
  
  obs <- compute_strata_probs(Y = Y, D = D, T = T, Z = Z)

      crit <- prepost_bounds_core(
        P = obs$P, V = obs$V, Q = obs$Q, W = obs$W,
        P_den = obs$P_den, V_den = obs$V_den,
        stable_mod = stable_mod,
        outcome_mono = outcome_mono,
        moderator_mono = moderator_mono,
        tau = tau,
        estimate_id_region = TRUE, N = N)

  no_con <- prepost_bounds_core(
        P = obs$P, V = obs$V, Q = obs$Q, W = obs$W,
        P_den = obs$P_den, V_den = obs$V_den,
        stable_mod = stable_mod,
        outcome_mono = outcome_mono,
        moderator_mono = moderator_mono,
        criterion_hat = crit$criterion, tau = tau, N = N)

  if (missing(g_at)) {
    rhos <- seq(0, no_con$max_rho, by = g_by)
  } else {
    rhos <- g_at
  }

  if (missing(t_at)) {
    thetas <- seq(0, no_con$max_theta, by = t_by)
  } else {
    thetas <- t_at
  }
  

  
  sens_out <- list()
  sens_out$rho <- rhos
  sens_out$theta <- thetas
  sens_out$lower <- matrix(NA, nrow = length(rhos), ncol = length(thetas))
  sens_out$upper <- matrix(NA, nrow = length(rhos), ncol = length(thetas))
  sens_out$criterion <- matrix(NA, nrow = length(rhos), ncol = length(thetas))
  ## sens_out$ms_p <- matrix(NA, nrow = length(rhos), ncol = length(thetas))
  sens_out$ci_lower <- matrix(NA, nrow = length(sens_out$rho),
                       ncol = length(sens_out$theta))
  sens_out$ci_upper <- matrix(NA, nrow = length(sens_out$rho),
                       ncol = length(sens_out$theta))

  for (r in seq_along(rhos)) {
    for (tt in seq_along(thetas)) {
      crit <- prepost_bounds_core(
        P = obs$P, V = obs$V, Q = obs$Q, W = obs$W,
        P_den = obs$P_den, V_den = obs$V_den,
        rho=rhos[r], theta = thetas[tt],
        stable_mod = stable_mod,
        outcome_mono = outcome_mono,
        moderator_mono = moderator_mono,
        tau = tau,
        estimate_id_region = TRUE, N = N)
      out <- prepost_bounds_core(
        P = obs$P, V = obs$V, Q = obs$Q, W = obs$W,
        P_den = obs$P_den, V_den = obs$V_den,
        rho= rhos[r], theta = thetas[tt],
        stable_mod = stable_mod,
        outcome_mono = outcome_mono,
        moderator_mono = moderator_mono,
        criterion_hat = crit$criterion, tau = tau, N = N)

      sens_out$lower[r, tt] <- out$lower
      sens_out$upper[r, tt] <- out$upper
      sens_out$criterion[r, tt] <- crit$criterion
      cns <- FALSE
      if (cns) {
        if (progress) cat("Finding upper CI via CNS...\n")
        upper_ci <- find_endpoint(
          out, obs, Y, D, T, Z, N, rhos[r],  thetas[tt],
          outcome_mono,
          moderator_mono, stable_mod, sims, upper = TRUE,
          conf_level, tau, progress)
      if (progress) cat("Finding lower CI via CNS...\n")
        lower_ci <- find_endpoint(
          out, obs, Y, D, T, Z, N, rhos[r],  thetas[tt],
          outcome_mono,
          moderator_mono, stable_mod, sims, upper = FALSE,
          conf_level, tau, progress)
        if (max(lower_ci$p_values) > 1 - conf_level) {
          sens_out$ci_lower[r, tt] <- min(setdiff(lower_ci$point_list, lower_ci$reject_list))
        } else {
          sens_out$ci_lower[r, tt] <- out$lower
        }
        if (max(upper_ci$p_values) > 1 - conf_level) {
          sens_out$ci_upper[r, tt] <- max(setdiff(upper_ci$point_list, upper_ci$reject_list))
        } else {
          sens_out$ci_upper[r, tt] <- out$upper
        }

      }
      

    }
  }
  
  ## ms_test <- array(NA, c(sims, length(sens_out$rho), length(sens_out$theta)))
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
        crit <- prepost_bounds_core(
          P = ostar$P, V = ostar$V, Q = ostar$Q, W = ostar$W,
          P_den = ostar$P_den, V_den = ostar$V_den,
          rho=rhos[r], theta = thetas[tt],
          stable_mod = stable_mod,
          outcome_mono = outcome_mono,
          moderator_mono = moderator_mono,
          estimate_id_region = TRUE, N = N)
        out <- prepost_bounds_core(
          P = ostar$P, V = ostar$V, Q = ostar$Q, W = ostar$W,
          P_den = obs$P_den, V_den = obs$V_den,
          rho= rhos[r], theta = thetas[tt],
          stable_mod = stable_mod,
          outcome_mono = outcome_mono,
          moderator_mono = moderator_mono,
          criterion_hat = crit$criterion, tau = 0.1, N = N)
        sens_lo[b, r, tt] <- out$lower
        sens_hi[b, r, tt] <- out$upper

        ## q_ms  <- prepost_bounds_inference(
        ##   P = obs$P, V = obs$V, Q = obs$Q,
        ##   P_star = ostar$P, V_star = ostar$V, Q_star = ostar$Q,
        ##   rho = rhos[r], theta = thetas[tt],
        ##   outcome_mono = outcome_mono,
        ##   moderator_mono = moderator_mono,
        ##   stable_mod = stable_mod,
        ##   criterion_hat = sens_out$criterion[r, tt], N = N, tau = 3)
        ## ms_test[b, r, tt] <- sens_out$criterion[r, tt] <= q_ms
      
      } 
    }
  }
  ## sens_out$ms_p <- apply(ms_test, c(2, 3), mean, na.rm = TRUE)
  sd_lo <- apply(sens_lo, c(2, 3), sd, na.rm = TRUE)
  sd_hi <- apply(sens_hi, c(2, 3), sd, na.rm = TRUE)
  for (i in seq_along(sens_out$rho)) {
    for (tt in seq_along(sens_out$theta)) {
      if (!is.na(sens_out$lower[i, tt])) {
        imbman <- imb.man.ci(sens_out$lower[i, tt], sens_out$upper[i, tt],
                             sd_lo[i, tt], sd_hi[i, tt], N = nrow(data),
                             alpha = conf_level)
        sens_out$ci_lower[i, tt] <- imbman[1]
        sens_out$ci_upper[i, tt] <- imbman[2]
      }
    }
  }

  out <- bind_cols(
    expand.grid(gamma = rhos, theta = thetas),
    lower = c(sens_out$lower),
    upper = c(sens_out$upper),
    ci_lower = c(sens_out$ci_lower),
    ci_upper = c(sens_out$ci_upper),
    criterion = c(sens_out$criterion)
    )
  return(out)
}




#' @importFrom Rglpk Rglpk_solve_LP
prepost_bounds_core <- function(P, V, Q, W, P_den, V_den,
                                moderator_mono = NULL,
                                outcome_mono = NULL, stable_mod = FALSE,
                                rho = 1, theta = 1,
                                estimate_id_region = FALSE,
                                criterion_hat = 0,
                                tau = 0.25, N,
                                fix = NULL) {

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
  slack_vars <- nrow(t_grid) + nrow(td_grid)
  num_strata <- nrow(psis)
  
  eq_params <- length(P) + length(V)
  ineq_params <- nrow(t_grid) + nrow(td_grid)
  num_params <- num_strata + 2 * eq_params + 2 * ineq_params + slack_vars    

  ## indices for extra variables
  P_mom <- seq_len(length(P)) + num_strata
  P_abs <- seq_len(length(P)) + max(P_mom)
  V_mom <- seq_len(length(V)) + max(P_abs)
  V_abs <- seq_len(length(V)) + max(V_mom)
  rho_mom <- seq_len(nrow(t_grid)) + max(V_abs)
  rho_abs <- seq_len(nrow(t_grid)) + max(rho_mom)
  rho_slack <- seq_len(nrow(t_grid)) + max(rho_abs) 
  theta_mom <- seq_len(nrow(td_grid)) + max(rho_slack)
  theta_abs <- seq_len(nrow(td_grid)) + max(theta_mom)
  theta_slack <- seq_len(nrow(td_grid)) + max(theta_abs)

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

  ## V_Y0_T_D0
  Veqmat <- matrix(0, nrow = length(V), ncol = num_params)
  V_levs <- list(Y0 = c(1, 0), T = c(1, 0), D0 = c(1, 0))
  V_vals <- expand.grid(V_levs)
  for (k in seq_along(V)) {
    pos <- which(psis$Y0 == V_vals$Y0[k] & psis$T == V_vals$T[k] &
                   psis$D0 == V_vals$D0[k])
    Veqmat[k, pos] <- 1
    Veqmat[k, V_mom[k]] <- 1
  }
  Vabsmat <- matrix(0, nrow = 2 * length(V), ncol = num_params)
  for (k in seq_along(P)) {
    Vabsmat[k, V_mom[k]] <- -1
    Vabsmat[k, V_abs[k]] <- 1
    Vabsmat[length(V) + k, V_mom[k]] <- 1
    Vabsmat[length(V) + k, V_abs[k]] <- 1
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
  
  ## theta constraints within TD0X strata
  theta_con <- matrix(0, nrow = nrow(td_grid), ncol = num_params)
  for (k in seq_len(nrow(td_grid))) {
    pos <- which(psis$Y1 != psis$Y0 & psis$T == td_grid$T[k] &
                   psis$D0 == td_grid$D0[k])
    theta_con[k, pos] <- 1
    theta_con[k, theta_slack[k]] <- 1
    theta_con[k, theta_mom[k]] <- 1
  }
  theta_abs_con <- matrix(0, nrow = 2 * nrow(td_grid), ncol = num_params)
  for (k in seq_len(nrow(t_grid))) {
    theta_abs_con[k, theta_mom[k]] <- -1
    theta_abs_con[k, theta_abs[k]] <- 1
    theta_abs_con[nrow(td_grid) + k, theta_mom[k]] <- 1
    theta_abs_con[nrow(td_grid) + k, theta_abs[k]] <- 1    
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
  zero_con <- matrix(0, nrow = length(zeros), ncol = num_params)
  for (k in seq_along(zeros)) {
    zero_con[k, zeros[k]] <- 1
  }
  zero_dir <- rep("==", times = length(zeros))
  zero_b <- rep(0, times = length(zeros))

  f_dir <- c(
    rep("==", nrow(Peqmat)),
    rep(">=", nrow(Pabsmat)),
    rep("==", nrow(Veqmat)),
    rep(">=", nrow(Vabsmat)),
    rep("==", nrow(sum_con)),
    zero_dir,
    rep("==", times = nrow(rho_con)),
    rep(">=", times = nrow(rho_abs_con)),
    rep("==", times = nrow(theta_con)),
    rep(">=", times = nrow(theta_abs_con))
  )

  b <- c(
    P,
    rep(0, times = nrow(Pabsmat)),
    V,
    rep(0, times = nrow(Vabsmat)),
    rep(1, times = nrow(sum_con)),
    zero_b,
    rep(rho, times = nrow(rho_con)),
    rep(0, times = nrow(rho_abs_con)),
    rep(theta, times = nrow(theta_con)),
    rep(0, times = nrow(theta_abs_con))
    )

  A <- rbind(
    Peqmat, Pabsmat, Veqmat, Vabsmat, sum_con,
    zero_con,
    rho_con, rho_abs_con, theta_con, theta_abs_con
  )

  unbounded <- c(P_mom, V_mom, rho_mom, theta_mom)
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
    f_obj[V_abs] <- sqrt(N)
    f_obj[rho_abs] <- sqrt(N)
    f_obj[theta_abs] <- sqrt(N)

    
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
    crit_con[V_abs] <- sqrt(N)
    crit_con[rho_abs] <- sqrt(N)
    crit_con[theta_abs] <- sqrt(N)
    A <- rbind(A, crit_con)
    b <- c(b, criterion_hat * (1 + tau))
    f_dir <- c(f_dir, "<=")
    if (all(is.finite(A)) & all(is.finite(b))) {
      ## u_res <- lpSolve::lp("max", f_obj, A, f_dir, b)
      ## l_res <- lpSolve::lp("min", f_obj, A, f_dir, b)
      ## lower <- ifelse(l_res$status == 0, l_res$objval, NA)
      ## upper <- ifelse(u_res$status == 0, u_res$objval, NA)
      
      ## ## grab each CATE here
      d1_pos <- which(psis$D0 == 1)
      d0_pos <- which(psis$D0 == 0)
      ## cate_u_1 <- sum(u_res$solution[d1_pos] * f_obj[d1_pos])
      ## cate_u_0 <- sum(u_res$solution[d0_pos] * f_obj[d0_pos])
      ## cate_l_1 <- sum(l_res$solution[d1_pos] * f_obj[d1_pos])
      ## cate_l_0 <- sum(l_res$solution[d0_pos] * f_obj[d0_pos])
      ## browser()
      u_res <- Rglpk::Rglpk_solve_LP(f_obj, A, f_dir, b, max = TRUE,
                                     bounds = lp_bounds,
                                     control = list(presolve = TRUE))
      l_res <- Rglpk::Rglpk_solve_LP(f_obj, A, f_dir, b, max = FALSE,
                                     bounds = lp_bounds,
                                     control = list(presolve = TRUE))
      lower <- ifelse(l_res$status == 0, l_res$optimum, NA)
      upper <- ifelse(u_res$status == 0, u_res$optimum, NA)
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

  }

  
  return(list(Q = Q, lower = lower, upper = upper,
              max_rho = max_rho, max_theta = max_theta,
              cate_1 = c(lower = cate_l_1, upper = cate_u_1),
              cate_0 = c(lower = cate_l_0, upper = cate_u_0),
              criterion = criterion, solution = solution,
              solution_L = solution_L, solution_U = solution_U))
}


prepost_bounds_inference <- function(P, V, Q, P_star, V_star, Q_star,
                                moderator_mono = NULL,
                                outcome_mono = NULL, stable_mod = FALSE,
                                rho = 1, theta = 1,
                                estimate_id_region = FALSE,
                                criterion_hat = 0,
                                tau = 0.25, N,
                                fix = NULL) {

  ## P is missing y0 and d0 so we add those
  ## P_Y1_D1_T
  P_levs <- list(Y1 = c(1, 0), D1 = c(1, 0), T = c(1, 0))
  P_vals <- expand.grid(P_levs)
  ## psi_Y1_Y0_D1_T_D0
  psi_levs <- list(Y1 = c(1, 0), Y0 = c(1, 0), D1 = c(1, 0), T = c(1, 0),
                   D0 = c(1, 0))
  psis <- expand.grid(psi_levs)

  P_diff <- P_star - P
  Q_diff <- Q_star - Q
  V_diff <- V_star - V

  ## deterministic constraints: sum to 1 in td_grid
  ## inequality contraints: gamma in t_grid, theta in td_grid  
  td_grid <- expand.grid(psi_levs[c("T", "D0")])
  t_grid <- expand.grid(T = c(1, 0))

  ## we need slack for both theta and gamma constraints
  slack_vars <- nrow(t_grid) + nrow(td_grid)
  num_strata <- nrow(psis)
  
  eq_params <- length(P) + length(V)
  ineq_params <- nrow(t_grid) + nrow(td_grid)
  num_params <- 2 * num_strata + 5 * eq_params + 5 * ineq_params + 2 * slack_vars    

  ## indices for extra variables
  H_par <- seq_len(num_strata) + num_strata
  P_mom <- seq_len(length(P)) + max(H_par)
  P_abs <- seq_len(length(P)) + max(P_mom)
  V_mom <- seq_len(length(V)) + max(P_abs)
  V_abs <- seq_len(length(V)) + max(V_mom)
  rho_mom <- seq_len(nrow(t_grid)) + max(V_abs)
  rho_abs <- seq_len(nrow(t_grid)) + max(rho_mom)
  rho_slack <- seq_len(nrow(t_grid)) + max(rho_abs) 
  theta_mom <- seq_len(nrow(td_grid)) + max(rho_slack)
  theta_abs <- seq_len(nrow(td_grid)) + max(theta_mom)
  theta_slack <- seq_len(nrow(td_grid)) + max(theta_abs)
  P_sam <- seq_len(length(P)) + max(theta_slack) 
  P_pos <- seq_len(length(P)) + max(P_sam)
  P_neg <- seq_len(length(P)) + max(P_pos)
  V_sam <- seq_len(length(V)) + max(P_neg) 
  V_pos <- seq_len(length(V)) + max(V_sam)
  V_neg <- seq_len(length(V)) + max(V_pos)
  rho_sam <- seq_len(nrow(t_grid)) + max(V_neg) 
  rho_pos <- seq_len(nrow(t_grid)) + max(rho_sam)
  rho_neg <- seq_len(nrow(t_grid)) + max(rho_pos)
  rho_h_slack <- seq_len(nrow(t_grid)) + max(rho_neg)
  theta_sam <- seq_len(nrow(td_grid)) + max(rho_h_slack) 
  theta_pos <- seq_len(nrow(td_grid)) + max(theta_sam)
  theta_neg <- seq_len(nrow(td_grid)) + max(theta_pos)
  theta_h_slack <- seq_len(nrow(td_grid)) + max(theta_neg)
  ## these will construct the moment conditions
  ## the {P/V}_ind will be the difference between the sum of the psi
  ## values and the corresponding P/V value. 
  Peqmat <- matrix(0, nrow = length(P), ncol = num_params)
  for (k in seq_along(P)) {
    pos1 <- which(psis$Y1 == P_vals$Y1[k] & psis$D1 == P_vals$D1[k] &
                    psis$T == P_vals$T[k] & psis$D0 == 1)
    pos0 <- which(psis$Y1 == P_vals$Y1[k] & psis$D1 == P_vals$D1[k] &
                    psis$T == P_vals$T[k] & psis$D0 == 0)
    Peqmat[k, pos1] <- sqrt(N) * Q_diff
    Peqmat[k, pos0] <- -sqrt(N) * Q_diff
    Peqmat[k, P_mom[k]] <- 1

    ## H variables
    Peqmat[k, pos1 + num_strata] <- Q
    Peqmat[k, pos0 + num_strata] <- 1 - Q
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
  ## these encode the absolute value transformation.
  ## they will have <= 0 so that the Pa_ind values
  ## will be greater than P_ind or -P_ind (ie, abs(P_ind))
  Psammat <- matrix(0, nrow = 2 * length(P), ncol = num_params)
  for (k in seq_along(P)) {
    pos1 <- which(psis$Y1 == P_vals$Y1[k] & psis$D1 == P_vals$D1[k] &
                    psis$T == P_vals$T[k] & psis$D0 == 1)
    pos0 <- which(psis$Y1 == P_vals$Y1[k] & psis$D1 == P_vals$D1[k] &
                    psis$T == P_vals$T[k] & psis$D0 == 0)

    Psammat[k, pos1] <- Q
    Psammat[k, pos0] <- 1 - Q
    Psammat[k, P_sam[k]] <- 1
    Psammat[length(P) + k, P_sam[k]] <- 1
    Psammat[length(P) + k, P_pos[k]] <- -1
    Psammat[length(P) + k, P_neg[k]] <- 1
  }
  P_sam_dir <- rep("==", times = nrow(Psammat))
  P_sam_b <- c(P, rep(0, times = length(P)))
  
  ## V_Y0_T_D0
  Veqmat <- matrix(0, nrow = length(V), ncol = num_params)
  V_levs <- list(Y0 = c(1, 0), T = c(1, 0), D0 = c(1, 0))
  V_vals <- expand.grid(V_levs)
  for (k in seq_along(V)) {
    pos <- which(psis$Y0 == V_vals$Y0[k] & psis$T == V_vals$T[k] &
                   psis$D0 == V_vals$D0[k])
    ## this drops out because it doesn't depend on the data
    ## Veqmat[k, pos] <- sqrt(N)
    Veqmat[k, V_mom[k]] <- 1
    Veqmat[k, pos + num_strata] <- 1
  }
  Vabsmat <- matrix(0, nrow = 2 * length(V), ncol = num_params)
  for (k in seq_along(P)) {
    Vabsmat[k, V_mom[k]] <- -1
    Vabsmat[k, V_abs[k]] <- 1
    Vabsmat[length(V) + k, V_mom[k]] <- 1
    Vabsmat[length(V) + k, V_abs[k]] <- 1
  }

  Vsammat <- matrix(0, nrow = 2 * length(V), ncol = num_params)
  for (k in seq_along(V)) {
    pos <- which(psis$Y0 == V_vals$Y0[k] & psis$T == V_vals$T[k] &
                   psis$D0 == V_vals$D0[k])
    Vsammat[k, pos] <- 1
    Vsammat[k, V_sam[k]] <- 1
    Vsammat[length(V) + k, V_sam[k]] <- 1
    Vsammat[length(V) + k, V_pos[k]] <- -1
    Vsammat[length(V) + k, V_neg[k]] <- 1
  }
  V_sam_dir <- rep("==", times = nrow(Vsammat))
  V_sam_b <- c(V, rep(0, times = length(V)))

  
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
    rho_con[k, pos10] <- -sqrt(N) * Q_diff
    rho_con[k, pos01] <- sqrt(N) * Q_diff
    rho_con[k, pos10 + num_strata] <- 1 - Q
    rho_con[k, pos01 + num_strata] <- Q

    ## slack variables
    rho_con[k, rho_h_slack[k]] <- 1
    rho_con[k, rho_mom[k]] <- 1
  }

  rho_abs_con <- matrix(0, nrow = 2 * nrow(t_grid), ncol = num_params)
  for (k in seq_len(nrow(t_grid))) {
    rho_abs_con[k, rho_mom[k]] <- -1
    rho_abs_con[k, rho_abs[k]] <- 1
    rho_abs_con[nrow(t_grid) + k, rho_mom[k]] <- 1
    rho_abs_con[nrow(t_grid) + k, rho_abs[k]] <- 1    
  }
  
  rho_sam_con <- matrix(0, nrow = 2 * nrow(t_grid), ncol = num_params)
  for (k in seq_len(nrow(t_grid))) {
    pos01 <- which(psis$D1 == 0 & psis$T == t_grid$T[k] & psis$D0 == 1)
    pos10 <- which(psis$D1 == 1 & psis$T == t_grid$T[k] & psis$D0 == 0)
    rho_sam_con[k, pos10] <- 1 - Q
    rho_sam_con[k, pos01] <- Q
    rho_sam_con[k, rho_sam[k]] <- 1
    rho_sam_con[k, rho_slack[k]] <- 1
    rho_sam_con[nrow(t_grid) + k, rho_sam[k]] <- 1
    rho_sam_con[nrow(t_grid) + k, rho_pos[k]] <- -1
    rho_sam_con[nrow(t_grid) + k, rho_neg[k]] <- 1
  }

  rho_sam_dir <- rep("==", times = 2 * nrow(t_grid))
  rho_sam_b <- c(rep(rho, times = nrow(t_grid)), rep(0, times = nrow(t_grid)))

  ## theta constraints within TD0X strata
  theta_con <- matrix(0, nrow = nrow(td_grid), ncol = num_params)
  for (k in seq_len(nrow(td_grid))) {
    pos <- which(psis$Y1 != psis$Y0 & psis$T == td_grid$T[k] &
                   psis$D0 == td_grid$D0[k])
    ## theta_con[k, pos] <- 1
    theta_con[k, pos + num_strata] <- 1
    theta_con[k, theta_h_slack[k]] <- 1
    theta_con[k, theta_mom[k]] <- 1
  }
  theta_abs_con <- matrix(0, nrow = 2 * nrow(td_grid), ncol = num_params)
  for (k in seq_len(nrow(t_grid))) {
    theta_abs_con[k, theta_mom[k]] <- -1
    theta_abs_con[k, theta_abs[k]] <- 1
    theta_abs_con[nrow(td_grid) + k, theta_mom[k]] <- 1
    theta_abs_con[nrow(td_grid) + k, theta_abs[k]] <- 1    
  }

  theta_sam_con <- matrix(0, nrow = 2 * nrow(td_grid), ncol = num_params)
  for (k in seq_len(nrow(td_grid))) {
    pos <- which(psis$Y1 != psis$Y0 & psis$T == td_grid$T[k] &
                   psis$D0 == td_grid$D0[k])
    theta_sam_con[k, pos] <- 1
    theta_sam_con[k, theta_sam[k]] <- 1
    theta_sam_con[k, theta_slack[k]] <- 1
    theta_sam_con[nrow(td_grid) + k, theta_sam[k]] <- 1
    theta_sam_con[nrow(td_grid) + k, theta_pos[k]] <- -1
    theta_sam_con[nrow(td_grid) + k, theta_neg[k]] <- 1
  }
  theta_sam_dir <- rep("==", times = 2 * nrow(td_grid))
  theta_sam_b <- c(rep(theta, times = nrow(td_grid)), rep(0, times = nrow(td_grid)))


  crit_con <- matrix(0, nrow = 1, ncol = num_params)
  crit_con[c(P_pos, P_neg)] <- sqrt(N)
  crit_con[c(V_pos, V_neg)] <- sqrt(N)
  crit_con[c(rho_pos, rho_neg)] <- sqrt(N)
  crit_con[c(theta_pos, theta_neg)] <- sqrt(N)
  crit_dir <- "<="
  crit_b <- criterion_hat * (1 + tau)
  
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
  zero_con <- matrix(0, nrow = length(zeros), ncol = num_params)
  for (k in seq_along(zeros)) {
    zero_con[k, zeros[k]] <- 1
    zero_con[k, H_par[zeros[k]]] <- 1/sqrt(N)
  }
  zero_dir <- rep("==", times = length(zeros))
  zero_b <- rep(0, times = length(zeros))

  pos_con <- cbind(
    diag(nrow = num_strata),
    matrix(0, nrow = num_strata, ncol = num_params - num_strata)
  ) 
  pos_dirs <- ifelse(seq_len(num_strata) %in% zeros, "==", ">=")
  pos_con_slack <- cbind(
    diag(nrow = num_strata),
    matrix(0, nrow = num_strata, ncol = num_params - num_strata)
  ) 
  pos_dirs <- ifelse(seq_len(num_strata) %in% zeros, "==", ">=")

  h_cons <- matrix(0, nrow = 2 * num_strata, ncol = num_params)
  for (k in seq_len(num_strata)) {
    h_cons[k, k] <- 1
    h_cons[k, H_par[k]] <- 1 / sqrt(N)
    h_cons[num_strata + k, k] <- 1
    h_cons[num_strata + k, H_par[k]] <- 1 / sqrt(N)
  }
  rho_slack_con <- matrix(0, nrow = nrow(t_grid), ncol = num_params)
  for (k in seq_len(nrow(t_grid))) {
    rho_slack_con[k, rho_slack[k]] <- 1
    rho_slack_con[k, rho_h_slack[k]] <- 1 / sqrt(N)
  }
  rho_slack_b <- rep(0, times = nrow(t_grid))
  rho_slack_dirs <- rep(">=", times = nrow(t_grid))
  
  theta_slack_con <- matrix(0, nrow = nrow(td_grid), ncol = num_params)
  for (k in seq_len(nrow(td_grid))) {
    theta_slack_con[k, theta_slack[k]] <- 1
    theta_slack_con[k, theta_h_slack[k]] <- 1 / sqrt(N)
  }
  theta_slack_b <- rep(0, times = nrow(td_grid))
  theta_slack_dirs <- rep(">=", times = nrow(td_grid))

  h_sum <- matrix(0, nrow = 1, ncol = num_params)
  h_sum[H_par] <- 1

  f_dir <- c(
    rep("==", nrow(Peqmat)),
    rep(">=", nrow(Pabsmat)),
    rep("==", nrow(Veqmat)),
    rep(">=", nrow(Vabsmat)),
    rep("==", nrow(sum_con)),
    zero_dir,
    rep("==", times = nrow(rho_con)),
    rep(">=", times = nrow(rho_abs_con)),
    rep("==", times = nrow(theta_con)),
    rep(">=", times = nrow(theta_abs_con)),
    rep("<=", times = num_strata),
    rep(">=", times = num_strata),
    P_sam_dir,
    V_sam_dir,
    rho_sam_dir,
    theta_sam_dir,
    crit_dir,
    rho_slack_dirs,
    theta_slack_dirs,
    "=="
  )

  b <- c(
    sqrt(N) * P_diff,
    rep(0, times = nrow(Pabsmat)),
    sqrt(N) * V_diff,
    rep(0, times = nrow(Vabsmat)),
    rep(1, times = nrow(sum_con)),
    zero_b,
    rep(rho, times = nrow(rho_con)),
    rep(0, times = nrow(rho_abs_con)),
    rep(theta, times = nrow(theta_con)),
    rep(0, times = nrow(theta_abs_con)),
    rep(1, times = num_strata),
    rep(0, times = num_strata),
    P_sam_b,
    V_sam_b,
    rho_sam_b,
    theta_sam_b,
    crit_b,
    rho_slack_b,
    theta_slack_b,
    0
    )

  A <- rbind(
    Peqmat, Pabsmat,
    Veqmat, Vabsmat,
    sum_con,
    zero_con,
    rho_con, rho_abs_con,
    theta_con, theta_abs_con,
    h_cons,
    Psammat,
    Vsammat, rho_sam_con, theta_sam_con,
    crit_con,
    rho_slack_con,
    theta_slack_con,
    h_sum
  )

  unbounded <- c(H_par, P_mom, V_mom, rho_mom, theta_mom,
                 P_sam, V_sam, rho_sam, theta_sam, rho_h_slack,
                 theta_h_slack)
  lp_bounds <- list(
    lower = list(ind = unbounded,
                 val = rep(-Inf, length(unbounded))),
    upper = list(ind = unbounded,
                 val = rep(Inf, length(unbounded)))
  )
  
  criterion <- NA
  if (!is.null(fix)) {
    fix_con <- matrix(0, nrow = 1, ncol = num_params)
    fix_con[which(psis$Y1 == 1 & psis$T == 1 & psis$D0 == 1)] <- 1
    fix_con[which(psis$Y1 == 1 & psis$T == 0 & psis$D0 == 1)] <- -1
    fix_con[which(psis$Y1 == 1 & psis$T == 1 & psis$D0 == 0)] <- -1
    fix_con[which(psis$Y1 == 1 & psis$T == 0 & psis$D0 == 0)] <- 1

    fix_h_con <- matrix(0, nrow = 1, ncol = num_params)
    fix_h_con[which(psis$Y1 == 1 & psis$T == 1 & psis$D0 == 1) + num_strata] <- 1
    fix_h_con[which(psis$Y1 == 1 & psis$T == 0 & psis$D0 == 1) + num_strata] <- -1
    fix_h_con[which(psis$Y1 == 1 & psis$T == 1 & psis$D0 == 0) + num_strata] <- -1
    fix_h_con[which(psis$Y1 == 1 & psis$T == 0 & psis$D0 == 0) + num_strata] <- 1

    A <- rbind(A, fix_con, fix_con, fix_h_con)
    f_dir <- c(f_dir, "<=", ">=", "==")
    b <- c(b, fix +0, fix - 0, 0)
  }

  f_obj <- rep(0, times = num_params)
  ## f_obj[c(P_pos, P_neg)] <- sqrt(N)
  ## f_obj[c(V_pos, V_neg)] <- sqrt(N)
  ## f_obj[c(rho_pos, rho_neg)] <- sqrt(N)
  ## f_obj[c(theta_pos, theta_neg)] <- sqrt(N)

  f_obj[P_abs] <- 1
  f_obj[V_abs] <- 1
  f_obj[rho_abs] <- 1
  f_obj[theta_abs] <- 1
  
  ## browser()
  min_crit <- Rglpk::Rglpk_solve_LP(
    f_obj, A, f_dir, b, max = FALSE,
    bounds = lp_bounds,
    control = list(presolve = TRUE, tm_limit = 10000)
    )
  ## if (is.na(criterion))  browser()
  criterion <- ifelse(min_crit$status == 0, min_crit$optimum, NA)
  ## AA <- rbind(Psammat,
  ##             sum_con,
  ##             fix_con)
  ## bb <- c(P_sam_b,
  ##         rep(1, times = nrow(sum_con)),
  ##         fix)
  ## dd <- c(P_sam_dir,
  ##         rep("==", nrow(sum_con)),
  ##         "==")
  ## min_crit <- Rglpk::Rglpk_solve_LP(
  ##   f_obj, AA, dd, bb, max = FALSE,
  ##   bounds = lp_bounds,
  ##   control = list(presolve = TRUE, verbose = TRUE)
  ##   )

  
  return(criterion)
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
