#' Data generation
#'
#' @description{Auxiliary function for generating simple artificial data sets with normally distributed coefficients
#' and regressors. Note that we only report this function for reproducibility of the simulations from the PhD thesis of the author.}
#'
#' @import mvtnorm
#' @import stats
#'
#' @param p Number of variables (columns).
#' @param n Number of observations (rows).
#' @param s Sparsity. Real number between 0 and 1. \code{s=1} (default) leads to a coefficient vector without zero entries.
#' @param xmean Mean of each of the normally distributed columns. Default is 0.
#' @param betamean Mean of each of the normally distributed coefficients. Default is 0.
#' @param betasd Standard deviation of the normally distributed coefficients. Default is 1.
#' @param snr Signal to noise ratio. Real number greater than zero. Default is 2.
#' @param rho Parameter for a Toeplitz covariance structure of the regressors. Real number between -1 and 1. Default is
#' 0 which corresponds to uncorrelated columns.
#'
#' @return \item{D}{Data matrix \eqn{(X,Y)}.}
#' \item{vars}{A list of the relevant variables.}
#' @export
#'
#' @examples genDataFromExamples(10,25,0.3)
genDataFromExamples<-function(p,n,s=1,xmean=0,betamean=0,betasd=1,snr=2,rho=0){

     if(missing(n)|missing(p)) stop("Number of rows or columns is missing.")
     if(p<=0) stop("p must be positive.")
     if(p!=round(p)) stop("p must be an integer.")
     if(n<=0) stop("n must be positive.")
     if(n!=round(n)) stop("n must be an integer.")
     if(s<=0|s>1) stop("s must be in ]0,1].")
     if(snr<=0) stop("Signal to noise ratio must be positive.")


     ## snr: Signal to noise ratio
     ## na: Number of randomly distributed NA's in the data matrix

     V<-diag(p)
     if(rho!=0){
             V<-matrix(numeric(p^2),ncol=p)
             V<-rho^(abs(col(V)-row(V)))
     }
     X<-rmvnorm(n,mean=rep(xmean,p),V)
     beta<-rnorm(p,betamean,betasd)
     ind<-1:p
     if(s!=1){
           ind<-sample(1:p,ceiling(s*p),replace=FALSE)
           beta[-ind]<-0
     }
     noise<-rnorm(n)
     Y<-X%*%beta
     k<-sqrt(var(Y)/(snr*var(noise)))
     Y<-Y+as.numeric(k)*noise
     D<-data.frame(X,Y)

     ## Output
     ## D: Data matrix (X,Y)
     ## vars: Position of the true relevant variables (i.e., with
     ## true non-zero coefficient)

     return(list("D"=D,"vars"=colnames(D)[sort(ind)]))
}
