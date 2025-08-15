#' Regression to the mean under bivariate \eqn{t}-distributions
#'
#' This function is based on the truncated bivariate \eqn{t}-distributions and creates confidence intervals and tests hypotheses under
#' adjustment of regression to the mean.
#'
#' The RTM effect is estimated from bivariate samples whose test statistics follows \eqn{t}-distributions truncated the baseline measurement
#' either \code{"left"} or \code{"right"}. The treatment effect is estimated and test under adjustment of RTM proposed by Umair et al. (2024).
#' The RTM expression for right cutoff and equal df is:
#' \deqn{R_r(c; n, \theta) =\frac{\sqrt{n} \, \Gamma\left(\frac{n + 1}{2}\right) \, \left(1 - \sin(\theta)\right)}{\sqrt{\pi} \, \Gamma\left(\frac{n}{2}\right)
#'  \, (n - 1) \, \left(1 + \frac{c^2}{n}\right)^{\frac{n-1}{2}} \, \left(1 - F_1(c)\right)}}
#' where \eqn{\Gamma(\cdot)} is the gamma function. Due to the symmetry of the \eqn{t}-distribution, the relation \eqn{1-F_1(c)=F_1(-c)} holds true.
#' For unequal degrees of freedom, a closed-form expression for the RTM does not exist, and the function must be evaluated numerically.
#'
#' @param data Numeric bivariate data in which the baseline measurement is above or below cutoff points.
#' @param df Vector giving the degree of freedoms of the variables
#' @param theta Numeric value of mixing angle
#' @param cutoff Numeric value shows from which the baseline measurement is above or below.
#' @param truncation character string, one of \code{"left"} or \code{"right"}
#' indicating the specification of the truncation. For left truncation distribution \code{"left"} refers to
#' that the true treatment effect and baseline measurements are below the cutoff point.
#' @param conf.level confidence level for the returned confidence interval,
#' restricted to lie between zero and one
#' @param tol Relative tolerance for numerical integration (smaller values give more precision but may increase computation time).
#' @param maxEval Maximum number of function evaluations for the numerical integrator (0 uses the integrator's default limit).
#' @param maxGrowIter Maximum number of iterations for the "growing-window" fallback integration when the main integration does not converge.
#' @param growFactor Multiplier controlling how quickly the integration window expands in the growing-window method.
#' @return A list of class \code{htest}, containing the following components:
#' \item{statistic}{the t-statistic, with names attribute \code{"t"}}
#' \item{p.value}{the p-value for the test}
#' \item{conf.int}{is a confidence
#' interval (vector of length 2) for the true  treatment mean. The
#' confidence level is recorded in the attribute \code{conf.level}}.
#' \item{estimate}{vector of length 2, giving the estimate of RTM and treatment effect.}
#' \item{null.value}{is the value of the treatment mean specified by the null hypothesis.}
#' \item{data.name}{a character string (vector of length 1) containing the actual names of the input vectors \code{data}}
#' @section Null Hypothesis: The null hypothesis is
#' that the difference between baseline measurement and followup measurement of the population from which \code{data} is drawn is zero.
#' For the one-sample t-tests, the null hypothesis is that the
#' population mean for difference of pre and post measurement in \code{data} is equal to zero.
#'
#' The alternative hypothesis in each case i.e \code{truncation="left"} or \code{truncation="right"} does not equal to zero.
#' @importFrom cubature adaptIntegrate
#' @importFrom stats t.test
#' @import PairedData
#' @export RTM.T
#' @author Muhammad Umair, Manzoor Khan
#' @seealso \code{\link{RTM.Norm}}
#' @references Shaw, W. T., & Lee, K. T. A. (2008). Bivariate Student t distributions with variable marginal degrees
#'  of freedom and independence. \emph{Journal of Multivariate Analysis}, 99(6), 1276-1287.
#'
#'  Umair, M., Khan, M., & Olivier, J. (2024). Accounting for regression to the mean under the bivariate t-distribution.
#'  \emph{Statistical Methods in Medical Research}, 33(9), 1624-1636. DOI:/10.1177/09622802241267808
#'
#' @examples
#' library(PairedData)
#' data("Anorexia")
#'   ### The z-score of cutoff point
#' cutoff<-(mean(Anorexia[,1])-95)/sd(Anorexia[,1])
#'   ### Testing of treatment means under RTM
#' theta<-0.67
#' RTM<-RTM.T(data=Anorexia,df=c(15,15),theta=theta,cutoff =  cutoff,truncation =  "left")
#' RTM
#'   ### Estimated RTM effect
#' RTM$estimate[1]
#'   ### Treatment effect (left cutoff), Total effect + RTM
#' RTM$estimate[2]
RTM.T <- function(data, df, theta, cutoff, truncation,
                  conf.level = 0.95,
                  tol = 1e-5,
                  maxEval = 0,
                  maxGrowIter = 120,
                  growFactor = 3) {

  ## -------------------------
  ## Argument checks
  ## -------------------------
  if (!is.matrix(data) && !is.data.frame(data))
    stop("`data` must be a matrix or data.frame with two columns (t1,t2).")
  if (ncol(data) < 2) stop("`data` must have at least two columns (t1, t2).")
  if (nrow(data) < 3) stop("The observations must be at least three.")
  if (length(df) != 2) stop("Provide exactly two degrees of freedom: df = c(df1, df2).")
  if (any(df <= 0)) stop("Degrees of freedom must be positive.")
  truncation <- match.arg(truncation, c("left", "right"))
  if (!is.numeric(conf.level) || length(conf.level) != 1 ||
      conf.level <= 0 || conf.level >= 1)
    stop("conf.level must be between 0 and 1 (exclusive).")

  ## -------------------------
  ## Main RTM calculation
  ## -------------------------

  if (df[1] == df[2]) {
    # Closed-form case
    n <- df[1]
    if (n <= 1) stop("Degrees of freedom must be > 1 for closed form.")

    x <- (n + cutoff^2) / n
    gamma_ratio <- exp(lgamma((n + 1)/2) - lgamma(n/2))
    denom_left <- if (truncation == "left") {
      stats::pt(cutoff, df = n)
    } else {
      1 - stats::pt(cutoff, df = n)
    }
    if (denom_left <= 0)
      stop("Truncation probability is zero or negative.")

    E.t1 <- (sqrt(n) * gamma_ratio) / (sqrt(pi) * denom_left * ((n - 1) * (x)^((n - 1)/2)))
    E.t2 <- E.t1 * sin(theta)

  } else {
    g1 <- function(t1, t2) {
      t1 * dbvt(t1, t2, df = df, theta = theta,
                cutoff = cutoff, truncation = truncation)
    }
    g2 <- function(t1, t2) {
      t2 * dbvt(t1, t2, df = df, theta = theta,
                cutoff = cutoff, truncation = truncation)
    }

    E.t1 <- compute_conditional_integral(
      g1, truncation, cutoff, df, tol, maxEval, maxGrowIter, growFactor
    )
    E.t2 <- compute_conditional_integral(
      g2, truncation, cutoff, df, tol, maxEval, maxGrowIter, growFactor
    )

    if (!is.finite(E.t1) || !is.finite(E.t2) || c(E.t1 && E.t2)==0)
      stop("Numerical integration failed to produce finite results.")

  }

  ## -------------------------
  ## Treatment test adjustment
  ## -------------------------
  t1 <- as.numeric(data[, 1])
  t2 <- as.numeric(data[, 2])

  if (truncation == "left") {
    RTM_val <- E.t1 - E.t2
    test_data <- t1 - t2 + RTM_val
  } else {
    RTM_val <- E.t1 - E.t2
    test_data <- t1 - t2 - RTM_val
  }

  RTM_val_rounded <- round(as.numeric(RTM_val), 4)
  test <- stats::t.test(test_data, conf.level = conf.level)

  estimates <- c(RTM = RTM_val_rounded,
                 Delta = round(as.numeric(test$estimate), 4))

  output <- list(
    statistic = test$statistic,
    parameter = test$parameter,
    p.value = test$p.value,
    alternative = "two.sided",
    estimate = estimates,
    null.value = c("treatment effect" = 0),
    conf.int = test$conf.int,
    method = "***TREATMENT MEANS UNDER REGRESSION TO THE MEAN***",
    data.name = sprintf("<%s: %d x %d>", class(data)[1], nrow(data), ncol(data))
  )

  class(output) <- "htest"
  output
}
