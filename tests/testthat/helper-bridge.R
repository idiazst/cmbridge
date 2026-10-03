make_small_bridge <- function(n = 20000, seed = 11) {
  set.seed(seed)
  B <- runif(n, -1, 1)
  V <- 0.97 * B + 0.03 * runif(n, -1, 1)
  beta <- 1.7 + 0.25 * V + 0.10 * V^2 - 0.08 * V^3
  M <- rbinom(n, 1, 1 / beta)
  Vobs <- V; Vobs[M == 0] <- NA_real_
  list(B = B, V = V, Vobs = Vobs, M = M)
}

truth_poly <- function(v) 1.7 + 0.25 * v + 0.10 * v^2 - 0.08 * v^3
