################################################
## Function to update strata probabilities, p ##
################################################

ps_tzd_crosswalk <- list()
ps_tzd_crosswalk[["111"]] <- c("s111", "s101", "s110", "s100")
ps_tzd_crosswalk[["011"]] <- c("s111", "s011", "s110", "s010")
ps_tzd_crosswalk[["101"]] <- c("s111", "s011", "s101", "s001")
ps_tzd_crosswalk[["001"]] <- c("s111", "s011", "s101", "s001")
ps_tzd_crosswalk[["110"]] <- c("s011", "s001", "s010", "s000")
ps_tzd_crosswalk[["010"]] <- c("s101", "s001", "s100", "s000")
ps_tzd_crosswalk[["100"]] <- c("s110", "s010", "s100", "s000")
ps_tzd_crosswalk[["000"]] <- c("s110", "s010", "s100", "s000")

update.p <- function(tr, z, d, omega, g, possible.strata) {

  ## create a PO for D strata grid and subset to possible strata
  ## this will help us map PO strata to observed strata
  strata_grid <- expand.grid(d11 = c(1, 0), d01 = c(1, 0), d0 = c(1, 0))
  strata_labs <- paste0("s", do.call(paste0, strata_grid))
  rownames(strata_grid) <- strata_labs
  strata_grid <- strata_grid[possible.strata, ]

  tzd_grid <- expand.grid(tr = c(1, 0), z = c(1, 0), d = c(1, 0))
  tzd_grid_str <- do.call(paste0, tzd_grid)
  tzd_data_str <- paste0(tr, z, d)
  n_levs <- nrow(tzd_grid)

  ## maps each TZD strata to the column of strata_grid corresponding
  ## to the PO of D observed in that strata. Allows us to find the
  ## "compatible strata" below
  ps_obs_crosswalk <- rep(NA, nrow(tzd_grid))
  ps_obs_crosswalk[tzd_grid$tr == 1 & tzd_grid$z == 1] <- 1
  ps_obs_crosswalk[tzd_grid$tr == 0 & tzd_grid$z == 1] <- 2
  ps_obs_crosswalk[tzd_grid$z == 0] <- 3
  
  p <- matrix(0, nrow = length(tr), ncol = length(possible.strata))
  colnames(p) <- possible.strata
  for (k in 1:n_levs) {
    k_units <- tzd_data_str == tzd_grid_str[k]
    ## which PO for D do we observe in this strata of TZD?
    which_d_po <- ps_obs_crosswalk[k]
    ## which PO strata are compatible with this level?
    
    compat_strata <- ps_tzd_crosswalk[[tzd_grid_str[k]]] ##which(strata_grid[, which_d_po] == tzd_grid$d[k])
    ## calculate numerators on the log scale
    compat_strata <- intersect(compat_strata, possible.strata)
    num <- omega[k_units, compat_strata, drop = FALSE] +
      g[k_units, compat_strata, drop = FALSE]
    ## logsumexp trick to avoid underflows
    ## see eg. https://blog.feedly.com/tricks-of-the-trade-logsumexp/
    num_max <- apply(num, 1, max)
    den <- log(rowSums(exp(num - num_max))) + num_max
    ## covert back to regular scale
    p[k_units, compat_strata] <- exp(as.matrix(num) - den)
    ## p[k_units, compat_strata] <- as.matrix(exp(num) / rowSums(exp(num)))
  }
  return(as.data.frame(p))
}


##################################################
## Function to update strata proportions, omega ##
##################################################

update.omega <- function(covars, psis) {
  num <- covars %*% psis
  num_max <- apply(num, 1, max)
  den <- log(rowSums(exp(num - num_max))) + num_max
  omega.mlogit <- num - den
  colnames(omega.mlogit) <- colnames(psis)
  omega.mlogit <- as.data.frame(omega.mlogit)
  return(omega.mlogit)
}

## strata should be a N x (J - 1) matrix of strata indicators with the
## omitted strata corresponding to the last column of current.psi
update.psi <- function(strata, X, m0, P0, current.psi) {
  N <- nrow(X)
  K <- ncol(X)
  J <- ncol(strata)
  n <- rep(1, N)
  psi <- current.psi
  kappa <- (strata - 0.5)
  for (j in seq_len(J)) {
    b0 <- P0[, , j] %*% m0[, j]
    Xb.notj <- X %*% psi[, -j]
    max.Xb <- apply(Xb.notj, 1, max)
    denom.j <- max.Xb + log(rowSums(exp(Xb.notj - max.Xb)))
    eta.j <- X %*% psi[, j] - denom.j

    w <- rpg.devroye(N, n, eta.j)
    P.like.j <- crossprod(X, X * w)
    b.like.j <- crossprod(X,  kappa[, j] + denom.j * w)

    P.post.j <- P.like.j + P0[, , j]
    V.post.j <- chol2inv(chol(P.post.j))
    m.post.j <- V.post.j %*% (b.like.j + b0)
    psi[, j] <- m.post.j + t(chol(V.post.j)) %*% rnorm(K)
  }
  psi
}


update.beta <- function(y, X, m0, P0, current.betas) {

  N <- nrow(X)
  K <- ncol(X)
  n <- rep(1, N)

  betas <- current.betas
  Xb <- as.numeric(X %*% betas)
  w <- BayesLogit::rpg.devroye(N, n, Xb)

  P.like <- crossprod(X, X * w)
  P.post <- P.like + P0
  S <- chol2inv(chol(P.post))
  b.post <- colSums(X * (y - 0.5)) + P0 %*% m0
  m.post <- S %*% b.post
  betas <- m.post + t(chol(S)) %*% rnorm(K)
  names(betas) <- colnames(X)
  betas
}
