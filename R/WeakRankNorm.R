#' Weak ranking family (normalized)
#'
#' @description{Gradient-free Gradient Boosting family for the normalized weak ranking loss function.}
#'
#' @import pcaPP
#' @param K Indicates that we are only interesting in the top \eqn{K} instances. Must be an integer between 1 and the number
#' \eqn{n} of observations.
#'
#' @details{A more intuitive loss function than the weak ranking loss thanks to its normalization to a maximum value
#' of 1. For example, if a number \eqn{c} of the top \eqn{K} instances has not been ranked at the top of the list, the
#' normalized weak ranking loss is \eqn{C/K}. \code{WeakRankNorm} returns a family object as in the package \code{mboost}.}
#'
#'
#' @return A Boosting family object
#' @export
#'
#' @references{Werner, T., Gradient-Free Gradient Boosting, PhD Thesis, Carl von Ossietzky University Oldenburg, 2020,
#' Remark (5.2.4)}
#' @references{T. Hothorn, P. Bühlmann, T. Kneib, M. Schmid, and B. Hofner. mboost: Model-Based
#' Boosting, 2017}
#'
WeakRankNorm<-function(K){
    if(missing(K)) stop("Need to specify the number K of best instances.")
    Family(
      ngradient=function(y,f) NULL,
      loss=function(y,f){
        if(missing(f)|missing(y)) stop("Need two vectors.")
        if(length(y)!=length(f)) stop("Objects lengths must be identical.")
        if(K!=round(K)) stop("K must be an integer.")
        if(K<=0) stop("K must be positive.")
        if(K>length(y)) stop("K must not be higher than the length of the vectors.")
        if(K==length(y)) warning("You are using K=n. The weak ranking loss will not produce a reasonable result.")
        indy<-order(y,decreasing=TRUE)
        indf<-order(f,decreasing=TRUE)
        s<-sum(indf[1:K]%in%indy[1:K])
        return((K-s)/K)
     },
     offset=function(y,f) NULL,
     name="Weak ranking family"
  )}
