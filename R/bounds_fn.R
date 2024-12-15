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







