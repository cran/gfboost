#' Weak ranking family
#'
#' @description{Gradient-free Gradient Boosting family for the weak ranking loss function.}
#'
#' @import pcaPP
#'
#' @param K Indicates that we are only interesting in the top \eqn{K} instances. Must be an integer between 1
#' and the number \eqn{n} of observations.
#'
#' @details{The weak ranking loss may be regarded as a classification loss. The parameter \code{K} defines the top of the list,
#' consisting of the best \eqn{K} instances according to their response values. Then the weak ranking loss penalizes
#' ''misclassification'' in the sense that instances belonging to the top of the list are ranked lower and vice versa.
#' \code{WeakRank} returns a family object as in the package \code{mboost}.}
#'
#'
#' @return  A Boosting family object
#' @export
#'
#' @references{Werner, T., Gradient-Free Gradient Boosting, PhD Thesis, Carl von Ossietzky University Oldenburg, 2020,
#' Remark (5.2.1)}
#' @references{T. Hothorn, P. Bühlmann, T. Kneib, M. Schmid, and B. Hofner. mboost: Model-Based
#' Boosting, 2017}
#'
#' @examples {y<-c(-3, 10.3,-8, 12, 14,-0.5, 29,-1.1,-5.7, 119)
#' yhat<-c(0.02, 0.6, 0.1, 0.47, 0.82, 0.04, 0.77, 0.09, 0.01, 0.79)
#' WeakRank(4)@risk(y,yhat)}
#' @examples {y<-c(-3, 10.3,-8, 12, 14,-0.5, 29,-1.1,-5.7, 119)
#' yhat<-c(0.02, 0.6, 0.1, 0.47, 0.82, 0.04, 0.77, 0.09, 0.01, 0.79)
#' WeakRank(5)@risk(y,yhat)}
WeakRank<-function(K){
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
       return(2*(K-s)/length(y))
       },
     offset=function(y,f) NULL,
     name="Weak ranking family"
    )}
