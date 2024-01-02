n <- 10
y <- rbinom(n, size = 1, prob = 0.5)
z <- rbinom(n, size = 1, prob = 0.5)
d <- rbinom(n, size = 1, prob = 0.5)
t <- rbinom(n, size = 1, prob = 0.5)
x1 <- rnorm(n)
x2 <- rnorm(n)
x3 <- rnorm(n)
df <- data.frame(y, z, d, t, x1, x2, x3)

test_that("prepost_gibbs runs", {
  out <- prepost_gibbs(y ~ t, data = df, prepost = ~ z, moderator = ~ d,
                       covariates = ~ x1 + x2 + x3, iter = 10)
  expect_type(out, "list")
  out_nc <- prepost_gibbs_nocovar(y ~ t, data = df, prepost = ~ z, moderator = ~ d,
                                  iter = 10)
  expect_type(out_nc, "list")
})
