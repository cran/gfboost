#' Hard ranking family
#'
#' @description{Gradient-free Gradient Boosting family for the hard ranking loss function including its fast computation.}
#'
#' @import pcaPP
#' @import mboost
#'
#' @details{The hard ranking loss is used to compare different orderings, usually the true ordering of instances of a data
#' set according to their responses with the predicted counterparts. The usage of the \code{pcaPP} package avoids the
#' cumbersome computation that would require \deqn{frac{n}{2(n-1)}} comparisons. \code{Rank} returns a family object
#' as in the package \code{mboost}.}
#' @return A Boosting family object
#' @export
#'
#' @references{Werner, T., Gradient-Free Gradient Boosting, PhD Thesis, Carl von Ossietzky University Oldenburg, 2020,
#' Equations (5.2.2) and (5.2.3)}
#' @references{T. Hothorn, P. Bühlmann, T. Kneib, M. Schmid, and B. Hofner. mboost: Model-Based
#' Boosting, 2017}
#'
#' @examples {y<-c(-3, 10.3,-8, 12, 14,-0.5, 29,-1.1,-5.7, 119)
#'  yhat<-c(0.02, 0.6, 0.1, 0.47, 0.82, 0.04, 0.77, 0.09, 0.01, 0.79)
#'  Rank()@risk(y,yhat)}
#' @examples {x<-1:6
#' z<-6:1
#' Rank()@risk(x,z)}
#' @examples {x<-1:6
#' z<-1:6
#' Rank()@risk(x,z)}
Rank<-function(){
     Family(
        ngradient=function(y,f) NULL,
        loss=function(y,f,w){
              if(missing(f)|missing(y)) stop("Need two vectors.")
              if(length(y)!=length(f)) stop("Objects lengths must be identical.")
              return((1-cor.fk(y,f))/2)
         },
        offset=function(y,f) NULL,
        name="Ranking family"
  )}
