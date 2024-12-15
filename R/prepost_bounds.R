## PREPOST BOUNDS ----

#' Run Prepost bounds
#'
#' @inheritParams post_bounds
#' @param prepost A one-sided formula with syntax ~ z, where z is the indicator
#'   variable for whether the moderator was measured pre- or post-treatment.
#' @param sims An integer specifying the number of bootstrap replications for
#'   the confidence intervals.
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
#' @param solver A character indicating what linear programming solver to use:
#'   "Rglpk" (the default) or "lpSolve".
#'
#' @return A list object containing bounds.
#'
#' @examples
#' data(land_experiment)
#' prepost_bounds(
#'   support ~ treat_comb,
#'   data = land_experiment,
#'   moderator = ~ land_insecure,
#'   prepost = ~ prepost,
#'   sims = 50
#' )
#' @export
prepost_bounds <- function(formula, data,  moderator,  prepost,
                           sims = 1000, conf_level = 0.95,
                           moderator_mono = NULL, outcome_mono = NULL,
                           stable_mod = FALSE,
                           tau = 0.25,
                           progress = TRUE, solver = "Rglpk") {

  outcome = all.vars(formula)[1]
  treat = all.vars(formula)[2]
  moderator = all.vars(moderator)[1]
  order_var = all.vars(prepost)[1]


  Y <- data[[outcome]]
  D <- data[[moderator]]
  T <- data[[treat]]
  Z <- data[[order_var]]
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
                               estimate_id_region = TRUE, N = N, solver = solver)
  bounds <- prepost_bounds_core(P = obs$P, V = obs$V, Q = obs$Q, W = obs$W,
                                P_den = obs$P_den, V_den = obs$V_den,
                                outcome_mono = om,
                               moderator_mono = mm,
                                stable_mod = stable_mod,
                               criterion_hat = crit$criterion, tau = tau, N = N, solver = solver)

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
      estimate_id_region = TRUE, N = N, solver = solver)
    bounds <- prepost_bounds_core(
      P = ostar$P, V = ostar$V, Q = ostar$Q, W = ostar$W,
      P_den = ostar$P_den, V_den = ostar$V_den,
      outcome_mono = om,
      moderator_mono = mm,
      stable_mod = stable_mod,
      criterion_hat = crit_star$criterion, tau = tau, N = N, solver = solver)
    ## q_ms <- prepost_bounds_inference(
    ##   P = obs$P, V = obs$V, Q = obs$Q,
    ##   P_star = ostar$P, V_star = ostar$V, Q_star = ostar$Q,
    ##   outcome_mono = om,
    ##   moderator_mono = mm,
    ##   stable_mod = stable_mod,
    ##   criterion_hat = crit$criterion, N = N,
    ##   tau= tau)
    ## ms_test[b] <- q_ms

    boot_hi[b] <- bounds$upper
    boot_lo[b] <- bounds$lower
  }
  if (empty_count > 0) {
    warning(empty_count, " bootstrap replications required resampling due to empty cells.", call. = FALSE)
  }
  ## out$ms_p <- mean(crit$criterion <= ms_test + 1e-6)
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





#' Run sensitivity analysis for the randomized moderator placement design
#'
#' @inheritParams post_sens
#' @param prepost A one-sided formula with syntax ~ z, where z is the indicator variable for whether the moderator was measured pre- or post-treatment. 
#' @param g_at Vector specifying what values to set the \eqn{\gamma} parameter
#'  to in the sensitivity analysis. Overrides \code{g_by}.
#' @param t_by Numeric indicating the grid spacing for the
#' \eqn{\theta} parameter that restricts what proportion of units have
#' their outcomes affected by the pre vs post-measurement of the
#' moderator.
#' @param t_at Vector specifying what values to set the \eqn{\theta} parameter
#'  to in the sensitivity analysis. Overrides \code{t_by}.
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
#' @param solver A character indicating what linear programming solver to use:
#'   "Rglpk" (the default) or "lpSolve".
#'
#' @return A list object containing sensitivity output.
#'
#' @examples
#' data(land_experiment)
##' prepost_sens(
##'   support ~ treat_comb,
##'   data = land_experiment,
##'   moderator = ~ land_insecure,
##'   prepost =  ~ prepost,
##'   g_by = 0.1,
##'   t_at = c(0.25, 1),
##'   sims = 50,
##'   moderator_mono = NULL
##' )
#' @export
prepost_sens <- function(formula, data, moderator, prepost,
                         g_by, g_at, t_by, t_at, sims = 1000, stable_mod = FALSE,
                         conf_level = 0.95, moderator_mono = NULL,
                         outcome_mono = NULL, tau = 0.25,
                         progress = TRUE, solver = "Rglpk") {


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
        estimate_id_region = TRUE, N = N, solver = solver)

  no_con <- prepost_bounds_core(
        P = obs$P, V = obs$V, Q = obs$Q, W = obs$W,
        P_den = obs$P_den, V_den = obs$V_den,
        stable_mod = stable_mod,
        outcome_mono = outcome_mono,
        moderator_mono = moderator_mono,
        criterion_hat = crit$criterion, tau = tau, N = N, solver = solver)

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
        estimate_id_region = TRUE, N = N, solver = solver)
      out <- prepost_bounds_core(
        P = obs$P, V = obs$V, Q = obs$Q, W = obs$W,
        P_den = obs$P_den, V_den = obs$V_den,
        rho= rhos[r], theta = thetas[tt],
        stable_mod = stable_mod,
        outcome_mono = outcome_mono,
        moderator_mono = moderator_mono,
        criterion_hat = crit$criterion, tau = tau, N = N, solver = solver)

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
          estimate_id_region = TRUE, N = N, solver = solver)
        out <- prepost_bounds_core(
          P = ostar$P, V = ostar$V, Q = ostar$Q, W = ostar$W,
          P_den = obs$P_den, V_den = obs$V_den,
          rho= rhos[r], theta = thetas[tt],
          stable_mod = stable_mod,
          outcome_mono = outcome_mono,
          moderator_mono = moderator_mono,
          criterion_hat = crit$criterion, tau = 0.1, N = N, solver = solver)
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

  out <- dplyr::bind_cols(
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
                                tau = 0.25, N, fix = NULL,
                                solver = "Rglpk") {

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
      
      ## ## grab each CATE here
      d1_pos <- which(psis$D0 == 1)
      d0_pos <- which(psis$D0 == 0)
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
