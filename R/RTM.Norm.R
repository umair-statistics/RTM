#' Regression to the mean under bivariate Normal distribution
#'
#' This function is based on the truncated bivariate Normal distributions and creates confidence intervals and tests hypotheses under
#' adjustment of regression to the mean
#'
#' The RTM effect is estimated from bivariate samples drawn from Normal distribution in which the baseline measurement
#' is truncated either \code{"left"} or \code{"right"}. The treatment effect is estimated and test under adjustment of RTM proposed by Khan and Oliver (2022). The expression
#' of RTM for right truncation is. \deqn{{R_r}\left( {{x_0};\Sigma } \right) = \left( {{\sigma _1} - \rho {\sigma _2}} \right)\frac{{\phi \left( {{z_0}} \right)}}{{1 - \Phi \left( {{z_0}} \right)}}}
#' Where \eqn{z_o} shows the z-score of cutoff point, \eqn{\phi(z_0)} is the density, and \eqn{\Phi(z_o)} is the CDF of standard normal distribution.
#'
#' @param data Numeric bivariate data in which the baseline measurement is above or below cutoff points
#' @param mu A vector giving the means of the variables.
#' @param Sigma a positive-definite symmetric matrix specifying the covariance matrix of the variables.
#' @param delta a single number representing the value of the difference in means specified by the null hypothesis
#' @param cutoff Numeric value shows from which the baseline measurement is above or below.
#' @param truncation character string, one of \code{"left"} or \code{"right"}
#' indicating the specification of the truncation. For left truncation distribution \code{"left"} refers to
#' that the true treatment effect and baseline measurements are below the cutoff point.
#' @param conf.level confidence level for the returned confidence interval,
#' restricted to lie between zero and one
#'
#' @return A list of class \code{htest}, containing the following components:
#' \item{statistic}{the t-statistic, with names attribute \code{"t"}}
#' \item{p.value}{the p-value for the test}
#' \item{conf.int}{is a confidence interval (vector of length 2) for the true  treatment mean. The
#' confidence level is recorded in the attribute \code{conf.level}.}
#' \item{estimate}{vector of length 2, giving the estimate of RTM and treatment effect.}
#' \item{null.value}{is the value of the treatment mean specified by the null hypothesis.}
#' \item{data.name}{a character string (vector of length 1) containing the actual names of the input vectors \code{data}}
#' @section Null Hypothesis: The null hypothesis is
#' that the difference between baseline measurement and followup measurement of the population from which \code{data} is drawn is equal to \code{delta}.
#' For the one-sample z-tests, the null hypothesis is that the
#' population mean for difference of pre and post measurement in \code{data} less that for \code{y} is \code{delta}.
#'
#' The alternative hypothesis is that the mean of difference between pre and post measurement is not equal to \code{delta}
#' @importFrom BSDA z.test
#' @importFrom stats dnorm pnorm
#' @author Manzoor and Oliver
#' @seealso \code{\link{RTM.T}}
#' @references Khan, M., & Olivier, J. (2019). Regression to the mean for the bivariate binomial distribution.
#' \emph{Statistics in medicine}, 38(13), 2391-2412.

#' @export RTM.Norm
#'
#' @examples
#' data("LeadExposedChildren")
#'   ### extract treatment group, baseline and followup measurement
#' treatGrp<- LeadExposedChildren[LeadExposedChildren$Group=="A",c(3,6)]
#'   ### cutoff point it is different for different studies
#' cutoff<-19.5
#'   ### define mean and variance covariance matrix
#' mu<-c(19.808,18.168)
#' varcov<-matrix(c(8.269^2,0.885*8.269*7.929,
#' 0.885*8.269*7.929,7.929^2),nrow = 2)
#'   ### estimate the RTM effect
#' RTM<-RTM.Norm(treatGrp,mu=mu,Sigma = varcov,cutoff = cutoff,truncation = "right")
#' RTM
#'   ### Total effect
#' mean(treatGrp[,1]-treatGrp[,2])
#'   ### treatment effect
#' RTM$estimate[2]
#'   ### RTM
#' RTM$estimate[1]

RTM.Norm <-
function(data, mu, Sigma, delta=0, cutoff, truncation,
    conf.level = 0.95)
{
  #Checking the conditions
  p <- length(mu)
  if (!all(dim(Sigma) == c(p, p)))
    stop("incompatible arguments")
  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values
  if (!all(ev >= -1e-06 * abs(ev[1L])))
    stop("'Sigma' is not positive definite")
  if (nrow(data)<3)
    stop("The observations must be at least three")
  choices <- c("left","right")
  alt <- pmatch(truncation, choices)
  truncation <- choices[alt]
  if (length(truncation) > 1 || is.na(truncation))
    stop("truncation must be one \"left\", \"right\"")
  if (!missing(conf.level))
    if (length(conf.level) != 1 || is.na(conf.level) ||
        conf.level < 0 || conf.level > 1)
      stop("conf.level must be a number between 0 and 1")
  if (truncation=="left") {
    mux<-mu[1]
    muy<-mu[2]
    sigx<-diag(Sigma)[1]
    sigy<-diag(Sigma)[2]
    cor<-(Sigma[2,1])/sqrt(sigx*sigy)
    z<-(mu[1]-cutoff)/sqrt(sigx)
    SigXY<-(dnorm(z)/pnorm(z))*(z-dnorm(z)/pnorm(z))*(sqrt(sigx)-
                                                        cor*sqrt(sigy))^2+sigx+sigy-2*cor*sqrt(sigx*sigy)
    RTM<-(sqrt(sigx)-sqrt(sigy)*cor)*dnorm(z)/pnorm(z)
    RTM<-round(RTM,4)
    #extract pre and post variables
    X<-data[,1]
    Y<-data[,2]
    #compute the treatment effect
    delta.est<-round(mean(X-Y+RTM),4)
    #test the treatment effect by eliminating RTM
    test<- z.test(X-Y+RTM,sigma.x =SigXY,mu=delta,conf.level=conf.level)
  }else{
    mux<-mu[1]
    muy<-mu[2]
    sigx<-diag(Sigma)[1]
    sigy<-diag(Sigma)[2]
    cor<-(Sigma[2,1])/sqrt(sigx*sigy)
    z<-(cutoff-mu[1])/sqrt(sigx)
    SigXY<-(dnorm(z)/(1-pnorm(z)))*(z-dnorm(z)/(1-pnorm(z)))*(sqrt(sigx)-
                                                                cor*sqrt(sigy))^2+sigx+sigy-2*cor*sqrt(sigx*sigy)

    RTM<-(sqrt(sigx)-sqrt(sigy)*cor)*dnorm(z)/(1-pnorm(z))
    RTM<-round(RTM,4)
    #extract pre and post variables
    X<-data[,1]
    Y<-data[,2]
    #compute the treatment effect
    delta.est<-round(mean(X-Y-RTM),4)
    #test the treatment effect by eliminating RTM
    test<- z.test(X-Y-RTM,sigma.x =SigXY,mu=delta,conf.level=conf.level)
  }

  #confidence interval
  cint<- test$conf.int
  #extract p-value, test statistics
  pval<- test$p.value
  z.stat<-test$statistic
  estimates<-c(RTM,delta.est)
  names(z.stat)<-"z"
  names(estimates)<-c("RTM","Delta")
  method<-"***TREATMENT MEANS UNDER REGRESSION TO THE MEAN***"
  mu<-delta
  names(mu)<-"treatment effect"
  dname <- deparse(substitute(data))
  output<-list(statistic = z.stat, p.value = pval,alternative="two.sided",
               estimate = estimates, null.value = mu, conf.int = cint,
               method = method,data.name = dname)

  attr(output, "class") <- "htest"
  return(output)
}
