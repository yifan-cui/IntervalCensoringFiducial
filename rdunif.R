rdunif <- function(n, b, a = 1) {
  stopifnot(is.numeric(a), length(a) == 1)
  stopifnot(is.numeric(b), length(b) == 1)
  
  a1 <- min(a, b)
  b1 <- max(a, b)
  
  sample(b1 - a1 + 1, n, replace = TRUE) + a1 - 1
}