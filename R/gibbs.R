
#' Run Gibbs sampler for the random moderator placement design
#'
#' @param formula A formula with syntax `y ~ t`, where `y` is the
#' (unquoted) name of the outcome and `t` is the (unquoted) name of the treatment.
#' @param data A data.frame containing variables in the formula, moderator, and covariates arguments.
#' @param prepost A one-sided formula with syntax `~ z`, where `z` is the
#' indicator variable for whether the moderator was measured  pre- or post-treatment.
#' @param moderator A one-sided formuala with syntax `~ d`, where `d`
#' is the (unquoted) name of the  moderator variable for the CATE.
#' @param covariates A one-sided formula with syntax `~ x1 + x2`,
#' where the right-hand side variables signify which covariates the
#' Gibbs will use to try and narrow the bounds.
#' @param iter Integer indicating the number of iterations for the
#' Gibbs sampler.
#' @param thin Integer indicating how often the Gibbs sampler should
#' save a draw.
#' @param burn Integer indicating how many iterations should be
#' performed before saving draws in the Gibbs sampler.
#' @param offset A numeric value indicating the center of the prior distribution for the covariate coefficients. 
#' @param monotonicity A logical signifying whether the model assumes monotonicity.
#' @param stable A logical signifying whether the model assumes that the
#' the pre vs post indicator does not affect the moderator under the
#' control condition for treatment.
#' @param saturated A logical indicating whether the coefficients on
#' the covariates are allowed to vary by the principal strata.
#' @param priors A list object containing the priors for the Gibbs
#' sampler. Priors include beta.precision, psi.precision, alpha,
#' y.alpha, and y.beta.
#'
#'
#' @return A list object containing Gibbs posterior quantities of interest and parameters.

#'
#' @examples
#' x <- "alfa,bravo,charlie,delta"

#' @importFrom progress progress_bar
#' @importFrom gtools rdirichlet
#' @importFrom BayesLogit rpg.devroye
#' @export
prepost_gibbs <- function(formula, data, prepost, moderator, covariates,
                          iter = 1000, thin = 1, burn = 0, offset = 0,
                          monotonicity = TRUE, stable = TRUE, saturated = TRUE,
                          priors) {


  # Y = binary outcome
  # T = treatment assignment
  # Z = post-test measurement of covariate indicator
  # D = binary covariate
  # covars = matrix of covariates
  # d1 = covariate under treatment + post-test
  # d0 = covariate under control + post-test
  # d* = covariate under pre-test
  # since we assume "stable moderator under control"
  # we must have d0 = d*
  # psi = strata proportion
  # n = strata number
  # mu = strata mean
  # g = strata conditional outcome distribution
  # iter = number of iterations
  # strata.interaction=T

  if (iter %% thin != 0) stop("iter must be multiple of thin")
  n.samp <- iter / thin

  if (missing(priors)) {
    priors <- generate.priors(covariates = TRUE)

  }

  if (prepost == ~1) {
    data$z <- rep(1, length(nrow(data)))
    prepost = ~z
  }

  ## TODO: add formula sanity checks
  lhs <- as.character(formula[2])
  if(missing(covariates)){
    all.forms <- list(formula, prepost, moderator)

  } else {
    all.forms <- list(formula, prepost, moderator, covariates)
  }
  big.rhs <- lapply(
    all.forms,
    function(x) as.character(x[length(x)])
  )


  big.rhs <- paste0(big.rhs, collapse = " + ")
  big.form <- as.formula(paste0(lhs, " ~ ", big.rhs))
  mf <- model.frame(big.form, data = data)

  ## TODO: add binary checks for tr/z/d
  ## NB: using tr and not t to avoid messing with transpose t()
  ## Here we use model.matrix to get the vectors instead of grabbing
  ## them from mf to avoid any issues like them being logicals
  ## TODO: make code robust to dropping intercept in these formulas
  tr <- model.matrix(formula, mf)
  tr.name <- colnames(tr)[ncol(tr)]
  tr <- tr[, ncol(tr)]
  z <- model.matrix(prepost, mf)
  z.name <- colnames(z)[ncol(z)]
  z <- z[, ncol(z)]
  d <- model.matrix(moderator, mf)[, 2]
  y <- model.response(mf)

  # take covariate formula entered as input and turn it into a model matrix
  # note: if there are any factor variables then model.matrix will expand them into dummies
  # covar.formula = as.formula(~ x1 + x2 + x3)
  if (missing(covariates)) {
    stop("no covariates specified. Use gibbs_nocovar() instead.")
  }
  
  covars <- model.matrix(covariates, mf)
  tz.mat <- cbind(tr, z, tr * z)
  colnames(tz.mat) <- c(tr.name, z.name, paste0(tr.name, ":", z.name))
  N <- nrow(mf)
  # save possible strata
  possible.strata <- c("s111", "s110", "s100", "s101", "s011", "s010",
                       "s001", "s000")

  mono.strata <- c("s111", "s110", "s100", "s010", "s000")
  stable.strata <- c("s111", "s100", "s011", "s000")

  if (monotonicity) {
    possible.strata <-  possible.strata[possible.strata %in% mono.strata]
  }
  if (stable) {
    possible.strata <- possible.strata[possible.strata %in% stable.strata]
  }

  ## output list
  out <- list()
  # how many betas do we need in the outcome model?
  # 1 for each covariate (with factors expanded into dummies)
  # 2 for each stratum, excluding one as the reference category,
  # for the dummy and the interaction with t
  # +3 for coefficients on t, z, and t:z
  if (!saturated) {
    betas.length <- dim(covars)[2] + 2 * (length(possible.strata) - 1) + 3
  } else {
    betas.length <- dim(covars)[2] + 4 * (length(possible.strata) - 1) + 3
  }

  out$betas <- matrix(NA, n.samp, betas.length)

  # Set initial values to betas vector, not in output matrix.
  # This allows us to more easily only every X iterations.
  betas <- rnorm(betas.length, mean = 0 + offset, sd = 5)

  # we'll use these names below as well
  strata.coef.names <- possible.strata[possible.strata != "s000"]
  if (!saturated) {
    colnames(out$betas) <- c(
      colnames(covars), tr.name, z.name, paste(tr.name, z.name, sep = ":"),
      strata.coef.names,
      paste0(tr.name, ":", strata.coef.names)
    )
  } else {
    colnames(out$betas) <- c(
      colnames(covars), tr.name, z.name, paste(tr.name, z.name, sep = ":"),
      strata.coef.names,
      paste0(tr.name, ":", strata.coef.names),
      paste0(z.name, ":", strata.coef.names),
      paste0(tr.name, ":", z.name, ':',strata.coef.names)
    )
  }

  names(betas) <- colnames(out$betas)
  betas_m0 <- rep(0, length(betas))
  betas_P0 <- diag(priors$beta.precision, nrow = length(betas))

  # create empty list of matrices to store the multinom coefs
  psis.length <- length(possible.strata)
  psis.width <- dim(covars)[2]
  psi.mat <- matrix(NA, n.samp, psis.width)
  colnames(psi.mat) <- colnames(covars)
  out$psis <- replicate(psis.length, psi.mat, simplify = FALSE)
  # add labels
  names(out$psis) <- possible.strata


  # Create a matrix of starting values for psis that conform to both
  # the mlogit function and the update omega function.
  total.psis <- dim(covars)[2] * psis.length
  psis.starting.vals <- rnorm(total.psis, mean = 0 + offset, sd = 5)
  psis <- matrix(psis.starting.vals, dim(covars)[2], psis.length)
  colnames(psis) <- possible.strata
  rownames(psis) <- colnames(covars)
  psis[, "s000"] <- 0


  ## output holder for the marginal proportion of the s strata
  out$s.tab <- matrix(0, nrow = n.samp, ncol = length(possible.strata))
  colnames(out$s.tab) <- possible.strata

  ## output holder for the deltas
  out$delta.1 <- rep(NA, n.samp)
  out$delta.2 <- rep(NA, n.samp)

  ## g = conditional outcome distribution
  g <- matrix(NA, nrow = N, ncol = length(possible.strata))
  colnames(g) <- possible.strata
  g <- as.data.frame(g)

  pb <- progress::progress_bar$new(total = iter + burn, clear = FALSE)
  pb$tick(0)

  for (i in 1:(iter + burn)) {
    pb$tick()
    ## cycle through all possible strata and calculate g (on log scale)
    for (s in possible.strata) {
      ## calculate log g at observed values of y, covars, t, z
      ## this side steps the need to calculate different g vectors
      ## for each t/z combo and then select them.
      this_mu <- get.mu(
        s, tr, z, covars, betas, possible.strata, type = "g")
      g[[s]] <- plogis(this_mu, log.p = TRUE) * y +
        plogis(this_mu, log.p = TRUE, lower.tail = FALSE) * (1 - y)
    }
    # current psis = multinomial coefficients
    # omegas = log strata proportions

    omega <- update_omega(covars, psis)

    ## update strata probabilities, based on updated g and current omega's
    ## NOTE: x1 and x2 do not explicitly enter this function
    ## but they do so implicitly because the omega's are selected based on x1, x2

    ## MB: this function now takes in g and omega and produces a
    ## matrix of the same size with columns corresponding to possible
    ## strata. Each row is the conditional posterior probability of
    ## each possible strata for that unit.
    p <- update_p(tr, z, d, omega, g, possible.strata)
    ## p$s000 <- 1 - rowSums(p[, colnames(p) != "s000"])

    ## sometimes there are floating point errors that we should round
    ## to zero.
    ## SH: i raised the multipler below to 50 instead of 2
    ## as i was still getting some small -ve values
    # p$s000 <- ifelse(abs(p$s000) < 2 * .Machine$double.eps, 0, p$s000)
    ## p$s000 <- ifelse(abs(p$s000) < 50 * .Machine$double.eps, 0, p$s000)

    ## if (any(p < 0)) browser()
    ## if (any(is.na(p))) browser()

    ## update imputed strata based on last probabilities
    strata <- apply(
      as.matrix(p),
      MARGIN = 1,
      FUN = function(pvec) {
        sample(x = possible.strata,
               size = 1,
               prob = pvec)
      }
    )

    ## strata <- sim_data$strata ## replace with true strata

    # if any covariates were entered as factors then they were expanded ino dummies by model.matrix
    # let's add those dummies (if any) into the main dataset
    ## covars.expanded.cols <- colnames(covars)[!colnames(covars) %in% colnames(data)]
    ## covars.expanded <- covars[, covars.expanded.cols]
    ## data.x <- cbind(data, covars.expanded)
    ## # we also need to "reduce" the data
    ## # instead of 1 row per observation
    ## # we want 1 row per unique combination of covariates in the regression
    ## # with an additional column tallying how many obs are in each combo
    ## strata.matrix <- data.x %>%
    ##   mutate(dummy=1) %>% # add column of 1's
    ##   group_by(strata) %>% mutate(temp_id = row_number()) %>% # add a temporary row_id (this is to workaround a silly error with "spread"...)
    ##   spread(key=strata, value=dummy, fill=0) %>% # create strata dummies
    ##   select(colnames(covars), starts_with("s")) %>%
    ##   group_by_at(colnames(covars)) %>% # group by covariate combinations
    ##   summarize_at(vars(starts_with("s")), mean, na.rm = TRUE) %>% # take strata means
    ##   ungroup()
    ## strata.matrix <- strata.matrix[,c(colnames(covars), possible.strata)]

    ## strata.outcome.matrix <- strata.matrix %>%
    ##   select(starts_with("s")) %>% select(-s000)
    ## strata.covar.matrix <- strata.matrix %>%
    ##   select(colnames(covars))
    ## nums <- data.x %>%
    ##   group_by_at(colnames(covars)) %>% # group by covariate combinations
    ##   tally()
    ## covar.nums <- nums$n
    ## create strata outcome matrix
    strata.mat <- create.indicators(strata, strata.coef.names)

    ## create the priors matrix (diagonal for each set set of
    ## coefficients)
    psis_P0 <- array(0, dim = c(ncol(covars), ncol(covars), ncol(strata.mat)))
    for (j in seq_len(ncol(strata.mat))) {
      psis_P0[, , j] <- diag(priors$psi.precision, ncol(covars))
    }

    psis <- update_psi(
      strata.mat, covars,
      m0 = array(0, dim = c(ncol(covars), ncol(strata.mat))),
      P0 = psis_P0,
      current.psi = psis
    )

    ## update the beta's (= binomial coefficients)
    ## create current strata imputation indicator/interaction matrices
    strata.int.mat <- tr * strata.mat
    colnames(strata.int.mat) <- paste0(tr.name, ":", colnames(strata.mat))



    if (saturated) {
      ## create z * strata matrix
      strata.int.mat.z <- z * strata.mat
      colnames(strata.int.mat.z) <- paste0(z.name, ":", colnames(strata.mat))

      ## create t * strata matrix
      strata.int.mat.z.t <- tr * z * strata.mat
      colnames(strata.int.mat.z.t) <- paste0(
        tr.name, ":", z.name, ":", colnames(strata.mat)
      )

      X <- cbind(
        covars,
        tz.mat,
        strata.mat,
        strata.int.mat,
        strata.int.mat.z,
        strata.int.mat.z.t
      )
    } else {
      X <- cbind(covars, tz.mat, strata.mat, strata.int.mat)
    }

    betas <- update_beta(
      y, X,
      m0 = betas_m0,
      P0 = betas_P0,
      current.betas = betas
    )


    ## DELTA METHOD 1 ----

    ## calculate "delta" for each iteration and store in output
    ## let's draw values of y11 and y01 for each observation based on the current betas
    ## but if the t and z values match for that observation, then we replace with the observed value y instead
    pred_y11 <- get.mu(strata, 1, 1, covars, betas, possible.strata)
    pred_y01 <- get.mu(strata, 0, 1, covars, betas, possible.strata)
    draw_y11 <- rbinom(N, 1, pred_y11)
    draw_y01 <- rbinom(N, 1, pred_y01)

    # replace with observed values where possible
    draw_y11[tr == 1 & z == 1] <- y[tr == 1 & z == 1]
    draw_y01[tr == 0 & z == 1] <- y[tr == 0 & z == 1]

    ## use list of strata corresponding to D*=1 and D*=0
    dstar.1.st <- c("s111", "s011", "s001", "s101")
    dstar.1.st <- dstar.1.st[dstar.1.st %in% possible.strata]
    dstar.0.st <- c("s110", "s010", "s000", "s100")
    dstar.0.st <- dstar.0.st[dstar.0.st %in% possible.strata]
    mu.11.dstar1 <- mean(draw_y11[strata %in% dstar.1.st])
    mu.01.dstar1 <- mean(draw_y01[strata %in% dstar.1.st])
    mu.11.dstar0 <- mean(draw_y11[strata %in% dstar.0.st])
    mu.01.dstar0 <- mean(draw_y01[strata %in% dstar.0.st])

    # print warnings
    if (sum(strata %in% dstar.1.st)==0){
      warning(paste("In iteration ", i, " there are no observations in strata", paste0(dstar.1.st, collapse=","), " so the in-sample delta is undefined."))
    }

    # print warnings
    if (sum(strata %in% dstar.0.st)==0){
      warning(paste("In iteration ", i, " there are no observations in strata", paste0(dstar.0.st, collapse=","), " so the in-sample delta is undefined."))
    }


    # calculate delta
    delta.1 <- (mu.11.dstar1 - mu.01.dstar1) - (mu.11.dstar0 - mu.01.dstar0)
    # if some strata are empty then delta.1 may end up being undefined
    # let's replace as a missing value
    delta.1 <- ifelse( is.nan(delta.1), NA, delta.1)

    ## DELTA METHOD 2 ----

    ## this is the "population estimate" where we compute delta for each observation using the predicted probabilities and then take an average

    ## convert betas to dataframe so we can pull out values based on coefficient name
    ## first get predicted mu's from the logit
    mu.11 <- vector("list", length = length(possible.strata))
    mu.01 <- vector("list", length = length(possible.strata))
    names(mu.11) <- possible.strata
    names(mu.01) <- possible.strata
    for (s in possible.strata) {
      mu.11[[s]] <- get.mu(s, 1, 1, covars, betas, possible.strata)
      mu.01[[s]] <- get.mu(s, 0, 1, covars, betas, possible.strata)
    }

    ## compute the relative strata probabilities
    ## no comma subsetting by dstar.*.st works b/c omega is df
    dstar1.prob <- rowSums(exp(omega[dstar.1.st]))
    rel.prob.dstar1 <- exp(omega[dstar.1.st]) / dstar1.prob
    dstar0.prob <- rowSums(exp(omega[dstar.0.st]))
    rel.prob.dstar0 <- exp(omega[dstar.0.st]) / dstar0.prob


    ## Average over the mus for those strata with the probs from above
    mu.11.dstar1 <- rowSums(mu.11[dstar.1.st] * rel.prob.dstar1)
    mu.01.dstar1 <- rowSums(mu.01[dstar.1.st] * rel.prob.dstar1)
    mu.11.dstar0 <- rowSums(mu.11[dstar.0.st] * rel.prob.dstar0)
    mu.01.dstar0 <- rowSums(mu.01[dstar.0.st] * rel.prob.dstar0)

    deltas <- (mu.11.dstar1 - mu.01.dstar1) - (mu.11.dstar0 - mu.01.dstar0)


    ## Save output if we're past burnin and on a thinning iteraction
    if ((i - burn) %% thin == 0 & i > burn) {
      ## Saved sample number/position
      samp <- (i - burn) / thin
      out$s.tab[samp, ] <- table(strata)[possible.strata] / N
      for (s in colnames(psis)) {
        out$psis[[s]][samp, ] <- psis[, s]
      }
      out$betas[samp, ] <- betas
      out$delta.1[samp] <- delta.1
      out$delta.2[samp] <- mean(deltas, na.rm = TRUE)

    }
  }

  return(out)

}



#' Run Gibbs sampler without covariates
#'
#' @param formula A formula with syntax `y ~ t`, where `y` is the name of the  outcome variable and `t` is the name of the treatment.
#' @param data A data frame containin the variables in the formula.
#' @param prepost A one-sided formula with syntax ~ z, where z is the indicator variable for whether the moderator was measured pre- or post-treatment. 
#' @param moderator A formuala with syntax ~ d, where d is the moderator variable for the CATE.
#' @param iter Numeric, number of iterations for the Gibbs
#' @param thin Numeric, thinning parameter for the Gibbs
#' @param burn Numeric, burn in rate for the Gibbs
#' @param monotonicity A logical signifying whether Gibbs assumes monotonicity.
#' @param stable A logical signifying whether Gibbs assumes stability.
#' @param priors A list object containing the priors for the Gibbs sampler. Priors include beta.precision, psi.precision, alpha, y.alpha, and y.beta.
#' @param predictive A logical indicator for whether to return prior predictive draws (`TRUE`) or posterior draws (`FALSE`, default).
#'
#' @return A list object containing Gibbs posterior quantities of interest and parameters.

#'
#' @examples
#' x <- "alfa,bravo,charlie,delta"

#' @importFrom progress progress_bar
#' @importFrom gtools rdirichlet
#' @importFrom BayesLogit rpg.devroye
#' @importFrom stats as.formula model.frame model.matrix model.response rnorm plogis rbinom rbeta dbeta ave runif
#' @export

prepost_gibbs_nocovar <- function(formula, data, prepost, moderator,
                                  iter = 1000, thin = 1, burn = 0,
                                  monotonicity = TRUE, stable = TRUE, priors,
                                  predictive = FALSE) {


  # Y = binary outcome
  # T = treatment assignment
  # Z = post-test measurement of covariate indicator
  # D = binary covariate
  # covars = matrix of covariates
  # d1 = covariate under treatment + post-test
  # d0 = covariate under control + post-test
  # d* = covariate under pre-test
  # since we assume "stable moderator under control"
  # we must have d0 = d*
  # psi = strata proportion
  # n = strata number
  # mu = strata mean
  # g = strata conditional outcome distribution
  # iter = number of iterations
  # strata.interaction=T

  if (iter %% thin != 0) stop("iter must be multiple of thin")
  n.samp <- iter / thin

  ## TODO: add formula sanity checks
  lhs <- as.character(formula[2])
  all.forms <- list(formula, prepost, moderator)
  big.rhs <- lapply(
    all.forms,
    function(x) as.character(x[length(x)])
  )

  use_data <- as.numeric(!predictive)


  big.rhs <- paste0(big.rhs, collapse = " + ")
  big.form <- as.formula(paste0(lhs, " ~ ", big.rhs))
  mf <- model.frame(big.form, data = data)

  ## TODO: add binary checks for tr/z/d
  ## NB: using tr and not t to avoid messing with transpose t()
  ## Here we use model.matrix to get the vectors instead of grabbing
  ## them from mf to avoid any issues like them being logicals
  ## TODO: make code robust to dropping intercept in these formulas
  tr <- model.matrix(formula, mf)
  tr.name <- colnames(tr)[ncol(tr)]
  tr <- tr[, ncol(tr)]
  z <- model.matrix(prepost, mf)
  z.name <- colnames(z)[ncol(z)]
  z <- z[, ncol(z)]
  d <- model.matrix(moderator, mf)[, 2]
  y <- model.response(mf)

  # take covariate formula entered as input and turn it into a model matrix
  # note: if there are any factor variables then model.matrix will expand them into dummies
  # covar.formula = as.formula(~ x1 + x2 + x3)
  tz.mat <- cbind(tr, z, tr * z)
  colnames(tz.mat) <- c(tr.name, z.name, paste0(tr.name, ":", z.name))

  N <- nrow(mf)
  # save possible strata
  possible.strata <- c("s111", "s110", "s100", "s101", "s011", "s010",
                       "s001", "s000")

  mono.strata <- c("s111", "s110", "s100", "s010", "s000")
  stable.strata <- c("s111", "s100", "s011", "s000")

  if (monotonicity) {
    possible.strata <-  possible.strata[possible.strata %in% mono.strata]
  }
  if (stable) {
    possible.strata <- possible.strata[possible.strata %in% stable.strata]
  }

  ## output list
  out <- list()
  # how many betas do we need in the outcome model?
  # 1 for each covariate (with factors expanded into dummies)
  # 2 for each stratum, excluding one as the reference category,
  # for the dummy and the interaction with t
  # +3 for coefficients on t, z, and t:z

  tzs_grid <- expand.grid(
    t = c(1, 0), z = c(1, 0),
    possible.strata = possible.strata
  )
  tzs_levs <- levels(interaction(tzs_grid, sep = "_"))
  mu <- runif(length(tzs_levs))
  names(mu) <- tzs_levs
  out$mu <- matrix(NA, n.samp, length(tzs_levs))
  colnames(out$mu) <- tzs_levs

  # create empty list of matrices to store the multinom coefs
  psi.length <- length(possible.strata)
  out$psis <- matrix(NA, n.samp, psi.length)
  # add labels
  colnames(out$psis) <- possible.strata

  # Create a matrix of starting values for psis that conform to both
  # the mlogit function and the update omega function.
  dstar.1.st <- c("s111", "s011", "s001", "s101")
  dstar.1.st <- dstar.1.st[dstar.1.st %in% possible.strata]
  dstar.0.st <- c("s110", "s010", "s000", "s100")
  dstar.0.st <- dstar.0.st[dstar.0.st %in% possible.strata]

  if (missing(priors)) {
    alpha <- rep(NA, psi.length)
    names(alpha) <- possible.strata
    alpha[dstar.1.st] <- 0.5 / length(dstar.1.st)
    alpha[dstar.0.st] <- 0.5 / length(dstar.0.st)
    priors <- list()
    priors$y.alpha <- 1 / (2 * psi.length)
    priors$y.beta <- 1 / (2 * psi.length)
  } else {
    alpha <- rep(priors$alpha, psi.length)
  }



  ## TODO: make these priors changeable
  ## if (monotonicity) alpha <- c(0.5, rep(0.5 / 4, length = 4))
  ## if (stable) alpha <- c(0.5, rep(0.5 / 3, length = 3))
  ## if (monotonicity & stable) alpha <- c(0.5, rep(0.5 / 2, length = 2))
  psis <- as.numeric(gtools::rdirichlet(1, alpha = alpha))
  names(psis) <- possible.strata

  ## output holder for the marginal proportion of the s strata
  out$s.tab <- matrix(0, nrow = n.samp, ncol = length(possible.strata))
  colnames(out$s.tab) <- possible.strata


  ## output holder for the deltas
  out$delta.1 <- rep(NA, n.samp)
  out$delta.2 <- rep(NA, n.samp)

  out$loglike <- rep(NA, n.samp)
  out$logpost <- rep(NA, n.samp)

  ## g = conditional outcome distribution
  g <- matrix(NA, nrow = N, ncol = length(possible.strata))
  colnames(g) <- possible.strata
  g <- as.data.frame(g)

  pb <- progress::progress_bar$new(total = iter + burn, clear = FALSE)
  pb$tick(0)
  for (i in 1:(iter + burn)) {
    pb$tick()
    ## cycle through all possible strata and calculate g (on log scale)
    for (s in possible.strata) {
      ## calculate log g at observed values of y, covars, t, z
      ## this side steps the need to calculate different g vectors
      ## for each t/z combo and then select them.

      this_mu <- paste0(tr, "_", z, "_", s)

      g[[s]] <- log(mu[this_mu]^y * (1-mu[this_mu]) ^ (1-y))## log(mu[this_mu]) * y +
      ## log(1 - mu[this_mu]) * (1 - y)
    }

    # current psis = multinomial coefficients
    # omegas = log strata proportions

    omega <- matrix(log(psis), nrow = N, ncol = psi.length, byrow = TRUE)
    colnames(omega) <- possible.strata
    ## update strata probabilities, based on updated g and current omega's
    ## NOTE: x1 and x2 do not explicitly enter this function
    ## but they do so implicitly because the omega's are selected based on x1, x2

    ## MB: this function now takes in g and omega and produces a
    ## matrix of the same size with columns corresponding to possible
    ## strata. Each row is the conditional posterior probability of
    ## each possible strata for that unit.
    if (predictive) {
      p <- psis
    } else {
      p <- update_p(tr, z, d, omega, g, possible.strata)
    }

    ## p$s000 <- 1 - rowSums(p[, colnames(p) != "s000"])

    ## sometimes there are floating point errors that we should round
    ## to zero.
    ## SH: i raised the multipler below to 50 instead of 2
    ## as i was still getting some small -ve values
    # p$s000 <- ifelse(abs(p$s000) < 2 * .Machine$double.eps, 0, p$s000)
    ## p$s000 <- ifelse(p$s000 < 0, 0, p$s000)

    ## if (any(is.na(p))) cat("psi: ", psis, "\nlog(psi): ", log(psis), "\n",
    ##                        "p: ", p, "min(g):", min(g), "\n")
    ## p[is.na(p)] <- 0

    ## if (any(is.na(p))) browser()

    ## update imputed strata based on last probabilities
    strata <- draw.strata(p, possible.strata)
    ## strata <- apply(
    ##   as.matrix(p),
    ##   MARGIN = 1,
    ##   FUN = function(pvec) {
    ##     sample(x = possible.strata,
    ##            size = 1,
    ##            prob = pvec)
    ##   }
    ## )

    ## strata <- sim_data$strata ## replace with true strata

    # if any covariates were entered as factors then they were expanded ino dummies by model.matrix
    # let's add those dummies (if any) into the main dataset
    ## covars.expanded.cols <- colnames(covars)[!colnames(covars) %in% colnames(data)]
    ## covars.expanded <- covars[, covars.expanded.cols]
    ## data.x <- cbind(data, covars.expanded)
    ## # we also need to "reduce" the data
    ## # instead of 1 row per observation
    ## # we want 1 row per unique combination of covariates in the regression
    ## # with an additional column tallying how many obs are in each combo
    ## strata.matrix <- data.x %>%
    ##   mutate(dummy=1) %>% # add column of 1's
    ##   group_by(strata) %>% mutate(temp_id = row_number()) %>% # add a temporary row_id (this is to workaround a silly error with "spread"...)
    ##   spread(key=strata, value=dummy, fill=0) %>% # create strata dummies
    ##   select(colnames(covars), starts_with("s")) %>%
    ##   group_by_at(colnames(covars)) %>% # group by covariate combinations
    ##   summarize_at(vars(starts_with("s")), mean, na.rm = TRUE) %>% # take strata means
    ##   ungroup()
    ## strata.matrix <- strata.matrix[,c(colnames(covars), possible.strata)]

    ## strata.outcome.matrix <- strata.matrix %>%
    ##   select(starts_with("s")) %>% select(-s000)
    ## strata.covar.matrix <- strata.matrix %>%
    ##   select(colnames(covars))
    ## nums <- data.x %>%
    ##   group_by_at(colnames(covars)) %>% # group by covariate combinations
    ##   tally()
    ## covar.nums <- nums$n

    ## create strata outcome matrix
    strata.mat <- create.indicators(strata, possible.strata)


    N.s <- colSums(strata.mat)
    psis <- gtools::rdirichlet(1, alpha = alpha + use_data * N.s)

    names(psis) <- possible.strata



    groups <- interaction(tr, z, strata, sep = "_")
    y_N1 <- y_N <- rep(0, length(tzs_levs))
    names(y_N1) <- names(y_N) <- tzs_levs
    y_N1[levels(groups)] <- tapply(y, groups, sum)
    y_N[levels(groups)] <- tapply(y, groups, length)
    y_N1[is.na(y_N1)] <- 0
    y_N[is.na(y_N)] <- 0
    mu <- rbeta(
      length(tzs_levs),
      shape1 = priors$y.alpha + use_data * y_N1,
      shape2 = priors$y.beta + use_data * (y_N - y_N1)
    )
    names(mu) <- tzs_levs

    psi_trans <- log(psis[1:(psi.length - 1)] / psis[psi.length])
    ll <- loglike_nocovar(
      c(psi_trans, log(mu / (1 - mu))),
      data = data.frame(t = tr, d = d, y = y, z = z),
      possible.strata = possible.strata
    )

    lprior <- log(gtools::ddirichlet(psis, alpha= alpha))
    lprior <- lprior + sum(dbeta(mu, priors$y.alpha, priors$y.beta, log = TRUE))
    lp <- -ll * use_data + lprior
    ## DELTA METHOD 1 ----

    ## This is the "in-sample" estimate where we use the current imputed strata to calculate various weighted means
    ## calculate "delta" for each iteration and store in output
    ## let's draw values of y11 and y01 for each observation based on the current betas
    ## but if the t and z values match for that observation, then we
    ## replace with the observed value y instead
    int_11 <- paste0("1_1_", strata)
    int_01 <- paste0("0_1_", strata)
    pred_y11 <- mu[int_11]
    pred_y01 <- mu[int_01]
    draw_y11 <- rbinom(N, 1, pred_y11)
    draw_y01 <- rbinom(N, 1, pred_y01)

    # replace with observed values where possible
    draw_y11[tr == 1 & z == 1] <- y[tr == 1 & z == 1]
    draw_y01[tr == 0 & z == 1] <- y[tr == 0 & z == 1]

    ## use list of strata corresponding to D*=1 and D*=0
    mu.11.dstar1 <- mean(draw_y11[strata %in% dstar.1.st])
    mu.01.dstar1 <- mean(draw_y01[strata %in% dstar.1.st])
    mu.11.dstar0 <- mean(draw_y11[strata %in% dstar.0.st])
    mu.01.dstar0 <- mean(draw_y01[strata %in% dstar.0.st])

    # print warnings
    if (sum(strata %in% dstar.1.st)==0){
      warning(paste("In iteration ", i, " there are no observations in strata", paste0(dstar.1.st, collapse=","), " so the in-sample delta is undefined."))
    }

    # print warnings
    if (sum(strata %in% dstar.0.st)==0){
      warning(paste("In iteration ", i, " there are no observations in strata", paste0(dstar.0.st, collapse=","), " so the in-sample delta is undefined."))
    }

    # calculate delta
    delta.1 <- (mu.11.dstar1 - mu.01.dstar1) - (mu.11.dstar0 - mu.01.dstar0)
    # if some strata are empty then delta.1 may end up being undefined
    # let's replace as a missing value
    delta.1 <- ifelse( is.nan(delta.1), NA, delta.1)

    ## DELTA METHOD 2 ----

    ## this is the "population estimate" where we compute delta for each observation using the predicted probabilities and then take an average

    ## convert betas to dataframe so we can pull out values based on coefficient name
    ## first get predicted mu's from the logit
    mu.11 <- rep(NA, length(possible.strata))
    mu.01 <- rep(NA, length(possible.strata))
    names(mu.11) <- possible.strata
    names(mu.01) <- possible.strata
    for (s in possible.strata) {
      mu.11[s] <- mu[paste0("1_1_", s)]
      mu.01[s] <- mu[paste0("0_1_", s)]
    }


    ## compute the relative strata probabilities
    rel.prob.dstar1 <- psis[dstar.1.st] / sum(psis[dstar.1.st])
    rel.prob.dstar0 <- psis[dstar.0.st] / sum(psis[dstar.0.st])

    ## Average over the mus for those strata with the probs from above
    mu.11.dstar1 <- sum(mu.11[dstar.1.st] * rel.prob.dstar1)
    mu.01.dstar1 <- sum(mu.01[dstar.1.st] * rel.prob.dstar1)
    mu.11.dstar0 <- sum(mu.11[dstar.0.st] * rel.prob.dstar0)
    mu.01.dstar0 <- sum(mu.01[dstar.0.st] * rel.prob.dstar0)

    deltas <- (mu.11.dstar1 - mu.01.dstar1) - (mu.11.dstar0 - mu.01.dstar0)


    ## Save output if we're past burnin and on a thinning iteraction
    if ((i - burn) %% thin == 0 & i > burn) {
      ## Saved sample number/position
      samp <- (i - burn) / thin
      s.tab <- table(strata)[possible.strata] / N
      s.tab[is.na(s.tab)] <- 0
      out$s.tab[samp, ] <- s.tab
      out$psis[samp, ] <- psis
      out$mu[samp, ] <- mu
      out$delta.1[samp] <- delta.1
      out$delta.2[samp] <- mean(deltas, na.rm=TRUE)
      out$loglike[samp] <- ll * use_data
      out$logpost[samp] <- lp
    }
  }

  return(out)

}




## Get predicted probability, mu ----


get.mu <- function(strata, t, z, covars, current.betas, possible.strata, type = "mu") {

  ## if given length 1 for any variable, make it compatible
  if (length(t) == 1) t <- rep(t, nrow(covars))
  if (length(z) == 1) z <- rep(z, nrow(covars))
  if (length(strata) == 1) strata <- rep(strata, nrow(covars))

  ## create strata indicator matrix
  used.strata <- possible.strata[possible.strata != "s000"]
  st.mat <- create.indicators(strata, used.strata)

  ## creat tz matrix
  tz.mat <- cbind(t = t, z = z, "t:z" = t * z)

  ## create t * strata matrix
  st.int.mat <- t * st.mat
  colnames(st.int.mat) <- paste0(colnames(st.mat), "*t")

  if('z:s111'%in%names(current.betas)){
    ## create z * strata matrix
    st.int.mat.z <- z * st.mat
    colnames(st.int.mat.z) <- paste0(colnames(st.mat), "*z")

    ## create t * strata matrix
    st.int.mat.z.t <- t * z * st.mat
    colnames(st.int.mat.z.t) <- paste0(colnames(st.mat), "*t*z")


  }

  ## combine into one design matrix
  if ("z:s111" %in% names(current.betas) | "s111*t*z" %in% names(current.betas)) {
    des.mat <- cbind(
      covars, tz.mat, st.mat,
      st.int.mat, st.int.mat.z, st.int.mat.z.t
    )
  } else{
    des.mat <- cbind(covars, tz.mat, st.mat, st.int.mat)
  }
  pred <- as.numeric(des.mat %*% current.betas)

  if (type == "mu") {
    mu.logit <- plogis(pred)
    mu.logit <- ifelse(is.na(mu.logit), 0, mu.logit)
    mu.logit <- unlist(mu.logit)
  } else if (type == "g") {
    mu.logit <- pred
  }
  return(mu.logit)
}



create.indicators <- function(x, vals) {
  ind.list <- lapply(vals, function(i) as.integer(x == i))
  out <- do.call(cbind, ind.list)
  colnames(out) <- vals
  return(out)
}

generate.priors <- function(covariates){
  out <- list()
  if (covariates) {
    out$beta.precision <- 0.1
    out$psi.precision <- 0.1
  } else {
    out$alpha <- 1
    out$y.alpha <- 1
    out$y.beta <- 1
  }
  out
}

draw.strata <- function(p, vals) {
  if (!is.matrix(p)) p <- as.matrix(p)
  u <- runif(nrow(p))
  c_p <- as.matrix(p) %*% upper.tri(diag(ncol(p)), diag = TRUE)
  inds_sel <- rowSums(u > c_p) + 1L
  vals[inds_sel]
}





loglike_nocovar <- function(theta, data, possible.strata, transformed = TRUE) {

  J <- length(possible.strata)
  N <- nrow(data)

  if (transformed) {
    psi <- exp(theta[1:(J - 1)]) / (1 + sum(exp(theta[1:(J - 1)])))
    psi <- c(psi, 1 / (1 + sum(exp(theta[1:(J - 1)]))))
    mu <- plogis(theta[J:length(theta)])
  } else {
    psi <- theta[1:J]
    mu <- theta[(J + 1):length(theta)]
  }
  tzs_grid <- expand.grid(
    t = c(1, 0), z = c(1, 0),
    possible.strata = possible.strata
  )
  tzs_levs <- levels(interaction(tzs_grid, sep = "_"))
  names(mu) <- tzs_levs

  g <- matrix(NA, nrow = N, ncol = J)
  colnames(g) <- possible.strata
  g <- as.data.frame(g)


  for (s in possible.strata) {
    ## calculate log g at observed values of y, covars, t, z
    ## this side steps the need to calculate different g vectors
    ## for each t/z combo and then select them.
    this_mu <- paste0(data$t, "_", data$z, "_", s)

    ## g[[s]] <- log(mu[this_mu]) * data$y + log(1 - mu[this_mu]) * (1 - data$y)
    g[[s]] <- ifelse(data$y == 1, mu[this_mu], 1 - mu[this_mu])
  }

  omega <- matrix(psi, nrow = N, ncol = J, byrow = TRUE)
  colnames(omega) <- possible.strata
  tzd_grid <- expand.grid(t = c(1, 0), z = c(1, 0), d = c(1, 0))
  tzd_grid_str <- do.call(paste0, tzd_grid)
  tzd_data_str <- with(data, paste0(t, z, d))
  n_levs <- nrow(tzd_grid)

  p <- matrix(0, nrow = N, ncol = length(possible.strata))
  colnames(p) <- possible.strata
  loglike <- 0
  for (k in 1:n_levs) {
    k_units <- tzd_data_str == tzd_grid_str[k]
    compat_strata <- ps_tzd_crosswalk[[tzd_grid_str[k]]]
    compat_strata <- intersect(compat_strata, possible.strata)

    num <- omega[k_units, compat_strata, drop = FALSE]  * g[k_units, compat_strata, drop = FALSE]
    num <- log(rowSums(num))
    loglike <- loglike + sum(num)

  }

  return(-loglike)
}
