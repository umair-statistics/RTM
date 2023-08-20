##' @name bv-t
##' @aliases rbvt
##' @aliases dbvt
##'
##' @title Bivariate \eqn{t}-distributions
##' @description
##' Generating random numbers, and density from bivariate \eqn{t}-distributions either truncated or complete.
##' The mean of complete distribution is equal to zero
##' @details
##' The bivariate \eqn{t}-distribution have identical marginals and equal or unequal degree of freedoms proposed by Shaw, W. T., & Lee, K. T. A. (2008).
##' The density function of bivariate \eqn{t}-distribution with equal degree of freedom is
##' \deqn{f\left( {t_1},{t_2};n,\theta \right) = \frac{1}{{2\pi \cos \left( \theta  \right)}}\frac{1}{{{{\left( {1 + \frac{\Delta }{n}} \right)}^{n/2 + 1}}}}}
##' where \eqn{-\infty < {t_1} < \infty}, \eqn{- \infty  < {t_2} < \infty}, \eqn{\Delta  = ({t_1^2 + t_2^2 - 2{t_1}{t_2}\sin (\theta )})^{-1}{{{\cos }^2}(\theta )}}, \eqn{n} is the degree of freedom, \eqn{\theta} is the mixing angle.
##'
##' @param t1 a numeric vector of baseline measurements
##' @param t2 a numeric vector of followup measurements
##' @param n the number of samples required
##' @param df a numeric vector giving the degree of freedoms of variables
##' @param theta a numeric value of mixing angle
##' @param cutoff a numeric value of cutoff point in which the baseline measurements are above or below
##' @param truncation a character string; one of \code{"left"} or \code{"right"} indicating the specification of distribution.
##' @return It give the matrix of random numbers and density whose dimension is determine by \code{n}, \eqn{t_1}, and \eqn{t_2}.
##' @importFrom stats pt rnorm rchisq
##' @importFrom hypergeo hypergeo
##' @importFrom graphics contour
##' @author Shaw, W. T., & Lee, K. T. A.
##' @seealso \code{\link{RTM.T}}
##' @references Shaw, W. T., & Lee, K. T. A. (2008). Bivariate Student t distributions
##'  with variable marginal degrees of freedom and independence. \emph{Journal of Multivariate Analysis}, 99(6), 1276-1287.
##'
##' @examples
##'   ### Generating a 1,000 samples
##' randomT<-rbvt(n=1000,df=c(15,9),theta=pi/4)
##' hist(randomT[,1])
##'
##' @rdname bv-t
##' @export
rbvt<-function (n,df,theta)
{
  #checking conditions
  if(length(df)!=2)
    stop("Proveide exactly two degree of freedoms")
  if(df[1]<=0||df[2]<=0)
    stop("Degree of freedom must be positive number")
  n1<-df[1]
  n2<-df[2]
  theta<-theta
  n<-n
  z1<-rnorm(n)
  z2<-rnorm(n)
  t1<-sqrt(n1/rchisq(n,n1))*z1
  t2<-sqrt(n2/rchisq(n,n2))*(z1*sin(theta)+z2*cos(theta))
  t<-cbind(t1,t2)
  return(round(t,9))
}
##'
##' @rdname bv-t
##' @examples
##'   #create bivariate t-distribution
##' t1    <- sort(randomT[,1])
##' t2    <- sort(randomT[,2])
##' df    <- c(15,9)
##' f     <- function(t1, t2) dbvt(t1,t2,df=c(15,9),theta= pi/4)
##' z     <- outer(t1, t2, f)
##'
##'   #create contour plot
##' contour(t1, t2, z)
##'
##'   ### Compute density from truncated distribution
##' randTrunc<-randomT[randomT[,1]>=0.5,]
##'   ### create right trucatated bivariate t-distribution
##' t1.Tr <- sort(randTrunc[,1])
##' t2.Tr <- sort(randTrunc[,2])
##' df    <- c(15,9)
##' f     <- function(t1.Tr, t2.Tr) dbvt(t1.Tr,t2.Tr,
##' df=c(15,9),theta= pi/4,cutoff = 0.5,truncation = "right")
##' z     <- outer(t1.Tr, t2.Tr, f)
##'
##'   #create contour plot
##' contour(t1.Tr, t2.Tr, z)
##' @export
dbvt <- function (t1,t2,df,theta,cutoff,truncation=NULL)
{
  #Checking conditions
  if (length(df)!=2)
    stop("Provide exactly two degrees of freedoms")
  if (df[1]<=0||df[2]<=0)
    stop("Degree of freedom must be positive number")
  if (!is.null(truncation)){
    choices <- c("left","right")
    alt <- pmatch(truncation, choices)
    truncation <- choices[alt]
    if (length(truncation) > 1|| is.na(truncation))
      stop("truncation must be one \"left\", \"right\"")
  }
  n1<-df[1]
  n2<-df[2]
  #compute the density
  alpha1<-1+t1^2/(n1*(cos(theta))^2)
  alpha2<-1+t2^2/(n2*(cos(theta))^2)
  gamma1<-(2*t1*t2*sin(theta))/(sqrt(n1*n2)*(cos(theta))^2)
  C<-1/(cos(theta)*pi*sqrt(n1*n2)*gamma(n1/2)*gamma(n2/2))

  a<-gamma((n1+1)/2)
  b<-gamma((n2+1)/2)
  ac<-(n1+1)/2
  bc<-(n2+1)/2
  x1<-gamma1^2/(4*alpha1*alpha2)
  c<-(n1/2+1)
  d<-(n2/2+1)
  dens<-C*alpha1^(-n1/2-1)*alpha2^(-n2/2-1)*(a*b*hypergeo(ac,bc,0.5,x1)*sqrt(alpha1*alpha2)+
                                               gamma1*gamma(n1/2+1)*gamma(n2/2+1)*hypergeo(n1/2+1,n2/2+1,3/2,x1))
  if (is.null(truncation)) {
    return(suppressWarnings(as.numeric(dens)))
  }
  else if (truncation == "left") {
    return(suppressWarnings(as.numeric(dens/pt(cutoff,n1))))
  }
  else {
    return(suppressWarnings(as.numeric(dens/(1-pt(cutoff,n1)))))
  }
}
