#' Localized ranking familty
#'
#' @description{Gradient-free Gradient Boosting family for the localized ranking loss function including its fast
#' computation.}
#' @import pcaPP
#' @import mboost
#' @param K Indicates that we are interesting in the top \eqn{K} instances and their correct ordering. Must be an integer
#' between 1 and the number \eqn{n} of observations.
#' @details{The localized ranking loss combines the hard and the weak ranking loss, i.e., it penalizes misrankings at
#' the top of the list (the best \eqn{K} instances according to the response value) and ''misclassification'' in the sense
#' that instances belonging to the top of the list are ranked lower and vice versa. The localized ranking loss already
#' returns a normalized loss that can take values between 0 and 1. \code{LocRank} returns a family object
#' as in the package \code{mboost}.}
#'
#' @return A Boosting family object
#' @export
#' @references{Werner, T., Gradient-Free Gradient Boosting, PhD Thesis, Carl von Ossietzky University Oldenburg, 2020,
#' Equation (5.2.5)}
#' @references{T. Hothorn, P. Bühlmann, T. Kneib, M. Schmid, and B. Hofner. mboost: Model-Based
#' Boosting, 2017}
#'
#' @examples {y<-c(-3, 10.3,-8, 12, 14,-0.5, 29,-1.1,-5.7, 119)
#'  yhat<-c(0.02, 0.6, 0.1, 0.47, 0.82, 0.04, 0.77, 0.09, 0.01, 0.79)
#'  LocRank(4)@risk(y,yhat)}
#' @examples {y<-c(-3, 10.3,-8, 12, 14,-0.5, 29,-1.1,-5.7, 119)
#'  yhat<-c(0.02, 0.6, 0.1, 0.47, 0.82, 0.04, 0.77, 0.09, 0.01, 0.79)
#'  LocRank(5)@risk(y,yhat)}
LocRank<-function(K){
  if(missing(K)) stop("Need to specify the number K of best instances.")
  Family(
    ngradient=function(y,f) NULL,
    loss=function(y,f,w){
      if(missing(f)|missing(y)) stop("Need two vectors.")
      if(length(y)!=length(f)) stop("Objects lengths must be identical.")
      if(K!=round(K)) stop("K must be an integer.")
      if(K<=0) stop("K must be positive.")
      if(K>length(y)) stop("K must not be higher than the length of the vectors.")
      if(K==length(y)) warning("You are using K=n, i.e., the hard ranking loss.")
      n<-length(y)
      yind<-y[rank(-f)<=K]
      find<-f[rank(-f)<=K]
      ##return(K*(K-1)/n/(n-1)*(1-cor.fk(yind,find))/2+(n-K)/n*WeakRank(K)@risk(y,f))  (unnormalized)
      return((K*(K-1)/n/(n-1)*(1-cor.fk(yind,find))/2+(n-K)/n*WeakRank(K)@risk(y,f))/(K*(K-1)/n/(n-1)+2*K*(n-K)/n^2))
    },
    offset=function(y,f) NULL,
    name="Localized ranking family"
  )}
