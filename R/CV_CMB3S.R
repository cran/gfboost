#' Cross-validated version of CMB-3S
#'
#' @description{Cross-validates the whole loss-based Stability Selection by aggregating several stable models according
#' to their performance on validation sets. Also computes a cross-validated test loss on a disjoint test set.}
#'
#' @param D Data matrix. Has to be an \eqn{n \times (p+1)-}dimensional data frame in the format \eqn{(X,Y)}. The \eqn{X-}part must not
#' contain an intercept column containing only ones since this column will be added automatically.
#' @param nsing Number of observations (rows) used for the SingBoost submodels.
#' @param Bsing Number of subsamples based on which the SingBoost models are validated. Default is 1. Not to confuse with parameter \code{B} for the Stability Selection.
#' @param B Number of subsamples based on which the CMB models are validated. Default is 100. Not to confuse with \code{Bsing} for CMB.
#' @param alpha Optional real number in \eqn{]0,1]}. Defines the fraction of best SingBoost models used in the aggregation step. Default is 1 (use all models).
#' @param singfam A SingBoost family. The SingBoost models are trained based on the corresponding loss function. Default is \code{Gaussian()} (squared loss).
#' @param evalfam A SingBoost family. The SingBoost models are validated according to the corresponding loss function. Default is \code{Gaussian()} (squared loss).
#' @param sing If \code{sing=FALSE} and the \code{singfam} family is a standard Boosting family that is contained in the package
#' \code{mboost}, the CMB aggregation procedure is executed for the corresponding standard Boosting models.
#' @param M An integer between 2 and \code{m_iter}. Indicates that in every \eqn{M-}th iteration, a singular iteration will be
#' performed. Default is 10.
#' @param m_iter Number of SingBoost iterations. Default is 100.
#' @param kap Learning rate (step size). Must be a real number in \eqn{]0,1]}. Default is 0.1 It is recommended to use
#' a value smaller than 0.5.
#' @param LS If a \code{singfamily} object that is already provided by \code{mboost} is used, the respective Boosting algorithm
#' will be performed in the singular iterations if \code{Ls} is set to \code{TRUE}. Default is \code{FALSE}.
#' @param best Needed in the case of localized ranking. The parameter \code{K} of the localized ranking loss will be
#' computed by \eqn{best \cdot n} (rounded to the next larger integer). Warning: If a parameter \code{K} is inserted into the
#' \code{LocRank} family, it will be ignored when executing SingBoost.
#' @param wagg Type of row weight aggregation. \code{'weights1'} indicates that the selection frequencies of the (best)
#' SingBoost models are averaged. \code{'weights2'} respects the validation losses for each model and downweights the ones
#' with higher validation losses.
#' @param gridtype Choose between \code{'pigrid'} and \code{'qgrid'}.
#' @param grid The grid for the thresholds (in \eqn{]0,1]}) or the numbers of final variables (positive integers).
#' @param ncmb Number of samples used for \code{CMB}. Integer that must be smaller than the number of samples in \code{D}.
#' @param CVind A list where each element contains a vector on length \eqn{n} (number of samples in the data matrix \code{D}) which
#' contains the strings \code{'tr'} (training set), \code{'v'} (validation set) and \code{'te'} (test set). This list can be easily
#' generated using the function \code{random_CVind}.
#' @param targetfam Target loss. Should be the same family as \code{evalfam}. Default is \code{Gaussian()} (squared loss).
#' @param print If set to \code{TRUE} (default), the number of the currently finished outer cross-validation loop is printed.
#' @param robagg Optional. If setting \code{robagg=TRUE}, the best SingBoost models are ignored when executing the
#' aggregation to avoid inlier effects. Only reasonable in combination with \code{lower}.
#' @param lower Optional argument. Only reasonable when setting \code{robagg=TRUE}. \code{lower} is a real number in \eqn{[0,1[} (a rather
#' small number is recommended) and indicates that the aggregation ignores the SingBoost models with the best
#' performances to avoid possible inlier effects.
#' @param singcoef Default is \code{FALSE}. Then the coefficients for the candidate stable models are computed by standard
#' linear regression (provided that the number of columns is smaller than the number of samples in the training set
#'  for each grid element). If set to \code{TRUE}, the coefficients are computed by SingBoost.
#' @param Mfinal Optional. Necessary if \code{singcoef=TRUE} to determine the frequency of singular iterations in the
#' SingBoost models.
#' @param ... Optional further arguments
#'
#' @details{In \code{CMB3S}, a validation set is given based on which the optimal stable model is chosen. The \code{CV.CMB3S}
#' function adds an outer cross-validation step such that both the training and the validation data sets (and
#' optionally the test data sets) are chosen randomly by disjointly dividing the initial data set. The aggregated
#' stable models form an ''ultra-stable'' model. It is strongly recommended to use this function is a parallelized
#' manner due to huge computation time.}
#'
#' @return \item{Cross-validated loss}{A vector containing the cross-validated test losses.}
#' \item{Ultra-stable column measure}{A vector containing the aggregated selection frequencies of the stable models.}
#' @export
#'
#' @references{Werner, T., Gradient-Free Gradient Boosting, PhD Thesis, Carl von Ossietzky University Oldenburg, 2020}
#'
CV.CMB3S<-function(D,nsing,Bsing=1,B=100,alpha=1,singfam=Gaussian(),evalfam=Gaussian(),sing=FALSE,M=10,m_iter=100,kap=0.1,LS=FALSE,best=1,wagg,gridtype,grid,ncmb,CVind,targetfam=Gaussian(),print=TRUE,robagg=FALSE,lower=0,singcoef=FALSE,Mfinal=10,...){

   if(missing(D)) stop("No training data argument.")
   if(missing(nsing)) stop("nsing is missing.")
   if(M!=round(M)) stop("M must be an integer.")
   if(M<=0) stop("M must be positive.")
   if(m_iter!=round(m_iter)) stop("m_iter must be an integer.")
   if(M>m_iter) stop("M must be at most as large as m_iter.")
   if(m_iter<=0) stop("m_iter must be positive.")
   if(kap<=0|kap>1) stop("Learning rate kap must be in ]0,1].")
   if(best<=0|best>1) stop("Parameter best must be in ]0,1].")
   if(nsing<=0) stop("nsing must be positive.")
   if(nsing!=round(nsing)) stop("nsing must be an integer.")
   if(nsing>dim(D)[1]) stop("nsing must not be higher than the number of rows in D.")
   if(Bsing<=0) stop("Bsing must be positive.")
   if(Bsing!=round(Bsing)) stop("Bsing must be an integer.")
   if(B<=0) stop("B must be positive.")
   if(B!=round(B)) stop("B must be an integer.")
   if(missing(grid)|missing(gridtype)) stop("grid or gridtype is missing.")
   if(missing(ncmb)) stop("ncmb is missing.")
   if(ncmb<=0) stop("ncmb must be positive.")
   if(ncmb!=round(ncmb)) stop("ncmb must be an integer.")
   if(ncmb>dim(D)[1]) stop("ncmb must not be higher than the number of rows in D.")
   if(ncmb<nsing) stop("ncmb must be at least as high as nsing.")
   if(Mfinal!=round(Mfinal)) stop("Mfinal must be an integer.")
   if(Mfinal<=0) stop("Mfinal must be positive.")

   loss<-numeric(length(CVind))
   ultrastab<-numeric(dim(D)[2])
   for(i in 1:length(CVind)){
         Dtrain<-D[which(CVind[[i]]=='tr'),]
         Dvalid<-D[which(CVind[[i]]=='v'),]
         Dtest<-D[which(CVind[[i]]=='te'),]
         res<-CMB3S(Dtrain=Dtrain,nsing=nsing,Bsing=Bsing,B=B,alpha=alpha,singfam=singfam,evalfam=evalfam,sing=sing,M=M,m_iter=m_iter,kap=kap,LS=LS,best=best,wagg=wagg,gridtype=gridtype,grid=grid,Dvalid=Dvalid,ncmb=ncmb,robagg=robagg,lower=lower,singcoef=singcoef,Mfinal=Mfinal)
         pred<-as.matrix(cbind(rep(1,dim(Dtest)[1]),Dtest[,res[[3]]]))%*%res[[1]]
         loss[i]<-targetfam@risk(pred,Dtest$Y)
         ultrastab<-(i-1)/i*ultrastab+1/i*res[[2]]
         if(print==TRUE) print(i)
   }

   return(list("Cross-validated loss"=loss,"Ultra-stable column measure"=ultrastab))
}
