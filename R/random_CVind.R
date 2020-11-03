#' Cross validation index generator
#'
#' @description{Simple auxiliary function for randomly generating the indices for training, validation and test data
#' for cross validation.}
#'
#' @param n Number of observations (rows).
#' @param ncmb Number of training samples for the SingBoost models in CMB. Must be an integer between 1 and \eqn{n}.
#' @param nval Number of validation samples in the CMB aggregation procedure. Must be an integer between 1 and \eqn{n-n_{cmb}-1}.
#' @param CV Number of cross validation steps. Must be a positive integer.
#'
#' @details{The data set consists of $n$ observations. \eqn{n_{cmb}} of them are used for the CMB aggregation procedure.
#' Note that within CMB itself, only a subset of these observations may be used for SingBoost training. The Stability
#' Selection is based on the validation set consisting of \eqn{n_{val}} observations. The cross-validated loss of the
#' final model is evaluated on the test data set with \eqn{n-n_{cmb}-n_{val}} observations. Clearly, all data sets need to
#' be disjoint.}
#'
#' @return \item{CVind}{List of row indices for training, validation and test data for each cross validation loop.}
#' @export
#'
random.CVind<-function(n,ncmb,nval,CV){
      if(n<=0) stop("n must be positive.")
      if(n!=round(n)) stop("n must be an integer.")
      if(ncmb<=0) stop("ncmb must be positive.")
      if(ncmb!=round(ncmb)) stop("ncmb must be an integer.")
      if(nval<=0) stop("nval must be positive.")
      if(nval!=round(nval)) stop("nval must be an integer.")
      if(CV<=0) stop("CV must be positive.")
      if(CV!=round(CV)) stop("CV must be an integer.")
      if(n-ncmb-nval<=0) stop("The sum of ncmb and nval is too large (has to be lower than n).")
      CVind<-list()
      for(k in 1:CV){
           ind<-numeric(n)
           indtr<-sample(1:n,ncmb)
           indval<-sample((1:n)[-indtr],nval)
           ind[indtr]<-'tr'
           ind[indval]<-'v'
           ind[-c(indtr,indval)]<-'te'
           CVind[[k]]<-ind
      }
      return(CVind)
}
