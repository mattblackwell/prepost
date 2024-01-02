set.seed(12345)
n <- 100
y <- rbinom(n, size = 1, prob = 0.5)
z <- rbinom(n, size = 1, prob = 0.5)
d <- rbinom(n, size = 1, prob = 0.5)
t <- rbinom(n, size = 1, prob = 0.5)
x1 <- rnorm(n)
x2 <- rnorm(n)
x3 <- rnorm(n)
df <- data.frame(y, z, d, t, x1, x2, x3)


test_that("prepost works", {
  np_prepost <- prepost_bounds(data = df,
                               formula = y ~ t,
                               moderator = ~ d,
                               prepost = ~ z,
                               sims = 10)
  expect_type(np_prepost, "list")
  np_post <- post_bounds(data = df[df$z == 1,],
                         formula = y ~ t,
                         moderator = ~ d,
                         sims = 10)
  expect_type(np_post, "list")
})
