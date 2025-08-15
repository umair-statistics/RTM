#' @noRd
safe_adapt_integrate <- function(f, lower, upper, tol_local, maxEval_local) {
  res <- tryCatch(
    cubature::adaptIntegrate(
      f = f,
      lowerLimit = lower,
      upperLimit = upper,
      tol = tol_local,
      maxEval = maxEval_local
    ),
    error = function(e) e
  )
  res
}

#' @noRd
# ---- helper: growing-window fallback integration ----
grow_window_integral <- function(f, truncation, cutoff,
                                 sd1, sd2,
                                 grow_tol,
                                 maxI,
                                 base_mult,
                                 tol_local,
                                 maxEval_local) {
  prev_val <- NA_real_
  B1 <- max(1, base_mult * sd1)
  B2 <- max(1, base_mult * sd2)

  for (I in seq_len(maxI)) {
    if (truncation == "left") {
      lower <- c(cutoff - I * B1, -I * B2)
      upper <- c(cutoff,           I * B2)
    } else {
      lower <- c(cutoff,           -I * B2)
      upper <- c(cutoff + I * B1,  I * B2)
    }

    res <- safe_adapt_integrate(f, lower, upper, tol_local, maxEval_local)
    if (inherits(res, "error")) next

    cur_val <- as.numeric(res$integral)

    if (is.na(prev_val)) {
      prev_val <- cur_val
      next
    }
    if (is.finite(cur_val) && is.finite(prev_val) &&
        abs(cur_val - prev_val) < grow_tol) {
      return(cur_val)
    }
    prev_val <- cur_val
  }
  NA_real_
}

#' @noRd
# ---- helper: compute conditional expectation via integration ----
compute_conditional_integral <- function(gfun, truncation, cutoff, df, tol, maxEval,
                                         maxGrowIter, growFactor) {
  integrand <- function(x) {
    val <- gfun(x[1], x[2])
    if (!is.finite(val)) val <- 0
    val
  }

  if (truncation == "left") {
    lower <- c(-Inf, -Inf)
    upper <- c(cutoff, Inf)
  } else {
    lower <- c(cutoff, -Inf)
    upper <- c(Inf, Inf)
  }

  res <- safe_adapt_integrate(integrand, lower, upper, tol, maxEval)
  if (!inherits(res, "error")) {
    val <- as.numeric(res$integral)
    if (is.finite(val) && !is.na(val)) return(val)
  }

  sd1 <- if (df[1] > 2) sqrt(df[1] / (df[1] - 2)) else 10
  sd2 <- if (df[2] > 2) sqrt(df[2] / (df[2] - 2)) else 10

  grow_window_integral(
    integrand, truncation, cutoff,
    sd1, sd2,
    grow_tol = max(tol, 1e-5),
    maxI = maxGrowIter,
    base_mult = growFactor,
    tol_local = tol,
    maxEval_local = maxEval
  )
}
