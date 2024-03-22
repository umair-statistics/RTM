#' Compute Regression to the mean under bivariate \eqn{t}-distributions
#'
#' This function is based on the truncated bivariate \eqn{t}-distributions and creates confidence intervals and tests hypotheses under
#' adjustment of regression to the mean.
#'
#' The RTM effect is estimated from bivariate samples drawn from \eqn{t}-distribution truncated the baseline measurement
#' either \code{"left"} or \code{"right"}. The treatment effect is estimated and test under adjustment of RTM.
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
#' @importFrom cubature pcubature
#' @importFrom stats t.test
#' @import PairedData
#' @export RTM.T
#' @author Muhammad Umair, Manzoor Khan
#' @seealso \code{\link{RTM.Norm}}
#' @references Shaw, W. T., & Lee, K. T. A. (2008). Bivariate Student t distributions with variable marginal degrees
#'  of freedom and independence. \emph{Journal of Multivariate Analysis}, 99(6), 1276-1287.
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
                  conf.level = 0.95)
{
  #Checking the conditions
  if (nrow(data)<3)
    stop("The observations must be at least three")
  if (length(df)!=2)
    stop("Provide exactly two degrees of freedoms")
  if (df[1]<0 || df[2]<0)
    stop("Degree of freedoms cannot be negative")
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
    if(df[1]==df[2]){
      #compute RTM from the derived formula
      n<-df[1]
      x<-(n+cutoff^(2))/n
      E.t1<-(sqrt(n)*gamma((n+1)/2))/(sqrt(pi)*gamma(n/2)*(pt(cutoff,n)))*(1/((n-1)*(x)^((n-1)/2)))
      E.t2<-(sqrt(n)*gamma((n+1)/2)*sin(theta))/(sqrt(pi)*gamma(n/2)*(pt(cutoff,n)))*(1/((n-1)*(x)^((n-1)/2)))
      RTM<-round(E.t1-E.t2,4)
      #Perform numerical integration if DoF are not equal
    }else{
      #E(t1|t1>c)
      RTM.t1 <- function(t,df,theta,cutoff,truncation){
        t[1]*dbvt(t[1],t[2],df,theta,cutoff,truncation)
      }
      #E(t2|t1>c)
      RTM.t2 <- function(t,df,theta,cutoff,truncation){
        t[2]*dbvt(t[1],t[2],df,theta,cutoff,truncation)
      }
      # Define the infinite limits
      lower_limit <- c(-Inf,-Inf)  # Lower limits for t1 and t2
      upper_limit <- c(cutoff, Inf)     # Upper limits for t1 and t2

      # Perform the double integration
      E.t1 <- pcubature(RTM.t1, lower_limit, upper_limit,df=df,theta=theta,cutoff=cutoff,truncation="left")$integral
      E.t2 <- pcubature(RTM.t2, lower_limit, upper_limit,df=df,theta=theta,cutoff=cutoff,truncation="left")$integral
      RTM<-round(E.t2-E.t1,4)
    }
    #extract pre and post variables
    t1<-as.matrix(data)[,1]
    t2<-as.matrix(data)[,2]
    #test the treatment effect by eliminating RTM
    test<- t.test(t1-t2+RTM,conf.level=conf.level)
    #compute the treatment effect
    delta<-round(test$estimate,4)
  }else{
    if(df[1]==df[2]){
      #compute RTM from the derived formula
      n<-df[1]
      x<-(n+cutoff^(2))/n
      E.t1<-(sqrt(n)*gamma((n+1)/2))/(sqrt(pi)*gamma(n/2)*(1-pt(cutoff,n)))*(1/((n-1)*(x)^((n-1)/2)))
      E.t2<-(sqrt(n)*gamma((n+1)/2)*sin(theta))/(sqrt(pi)*gamma(n/2)*(1-pt(cutoff,n)))*(1/((n-1)*(x)^((n-1)/2)))
      RTM<-round(E.t1-E.t2,4)
    }else{
      #E(t1|t1>c)
      RTM.t1 <- function(t,df,theta,cutoff,truncation){
        t[1]*dbvt(t[1],t[2],df,theta,cutoff,truncation)
      }
      #E(t2|t1>c)
      RTM.t2 <- function(t,df,theta,cutoff,truncation){
        t[2]*dbvt(t[1],t[2],df,theta,cutoff,truncation)
      }
      # Define the infinite limits
      lower_limit <- c(cutoff,-Inf)  # Lower limits for t1 and t2
      upper_limit <- c(Inf, Inf)     # Upper limits for t1 and t2

      # Perform the double integration
      E.t1 <- pcubature(RTM.t1, lower_limit, upper_limit,df=df,theta=theta,cutoff=cutoff,truncation="right")$integral
      E.t2 <- pcubature(RTM.t2, lower_limit, upper_limit,df=df,theta=theta,cutoff=cutoff,truncation="right")$integral
      RTM<-round(E.t1-E.t2,4)
    }
    #extract pre and post variables
    t1<-data[,1]
    t2<-data[,2]
    #test the treatment effect by eliminating RTM
    test<- t.test(t1-t2-RTM,conf.level=conf.level)
    #compute the treatment effect
    delta<-round(test$estimate,4)
  }

  #confidence interval
  cint<- test$conf.int
  #extract p-value, test statistics
  pval<- test$p.value
  t.stat<-test$statistic
  df<-test$parameter
  estimates<-c(RTM,delta)
  names(t.stat)<-"t"
  names(estimates)<-c("RTM","Delta")
  method<-"***TREATMENT MEANS UNDER REGRESSION TO THE MEAN***"
  mu<-0
  names(mu)<-"treatment effect"
  dname <- deparse(substitute(data))
  output<-list(statistic = t.stat,parameter=df, p.value = pval,alternative="two.sided",
               estimate = estimates, null.value = mu, conf.int = cint,
               method = method,data.name = dname)

  attr(output, "class") <- "htest"
  return(output)
}
