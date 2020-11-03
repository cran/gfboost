#' CMB stability paths
#'
#' @description{Draws a Stability plot for CMB.}
#'
#' @import stats
#' @import graphics
#'
#' @param D Data matrix. Has to be an \eqn{n \times (p+1)-}dimensional data frame in the format \eqn{(X,Y)}. The \eqn{X-}part must not
#' contain an intercept column containing only ones since this column will be added automatically.
#' @param nsing Number of observations (rows) used for the SingBoost submodels.
#' @param Bsing Number of subsamples based on which the SingBoost models are validated. Default is 1. Not to confuse with parameter \code{B} for the Stability Selection.
#' @param alpha Optional real number in \eqn{]0,1]}. Defines the fraction of best SingBoost models used in the aggregation step. Default is 1 (use all models).
#' @param singfam A SingBoost family. The SingBoost models are trained based on the corresponding loss function. Default is \code{Gaussian()} (squared loss).
#' @param evalfam A SingBoost family. The SingBoost models are validated according to the corresponding loss function. Default is \code{Gaussian()} (squared loss).
#' @param sing If \code{sing=FALSE} and the \code{singfam} family is a standard Boosting family that is contained in the package
#' \code{mboost}, the CMB aggregation procedure is executed for the corresponding standard Boosting models.
#' @param Mseq A vector of different values for \eqn{M}.
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
#' @param robagg Optional. If setting \code{robagg=TRUE}, the best SingBoost models are ignored when executing the
#' aggregation to avoid inlier effects. Only reasonable in combination with \code{lower}.
#' @param lower Optional argument. Only reasonable when setting \code{robagg=TRUE}. \code{lower} is a real number in \eqn{[0,1[} (a rather
#' small number is recommended) and indicates that the aggregation ignores the SingBoost models with the best
#' performances to avoid possible inlier effects.
#' @param B Number of subsamples of size \eqn{n_{cmb}} of the training data for CMB aggregation.
#' @param ncmb Number of samples used for \code{CMB}. Integer that must be smaller than the number of samples in \code{Dtrain}.
#' @param ... Optional further arguments
#'
#' @return \item{relev}{List of relevant variables (represented as their column number).}
#' \item{ind}{Vector of relevant variables (represented as their column number).}
#' @export
#'
#' @references{Werner, T., Gradient-Free Gradient Boosting, PhD Thesis, Carl von Ossietzky University Oldenburg, 2020}
#'
CMB.stabpath<-function(D,nsing,Bsing=1,alpha=1,singfam=Gaussian(),evalfam=Gaussian(),sing=FALSE,Mseq,m_iter=100,kap=0.1,LS=FALSE,best=1,wagg,robagg=FALSE,lower=0,B,ncmb,...){
        if(missing(D)) stop("No data argument.")
        if(missing(nsing)|missing(ind)) stop("nsing or ind is missing.")
        if(any(Mseq!=round(Mseq))) stop("Mseq must only contain integers.")
        if(any(Mseq<=0)) stop("Mseq must only contain positive values.")
        if(m_iter!=round(m_iter)) stop("m_iter must be an integer.")
        if(any(Mseq>m_iter)) stop("Mseq must only contain values that are at most as large as m_iter.")
        if(m_iter<=0) stop("m_iter must be positive.")
        if(kap<=0|kap>1) stop("Learning rate kap must be in ]0,1].")
        if(best<=0|best>1) stop("Parameter best must be in ]0,1].")
        if(nsing<=0) stop("nsing must be positive.")
        if(nsing!=round(nsing)) stop("nsing must be an integer.")
        if(nsing>=dim(D)[1]) stop("nsing must be smaller than the number of rows in D.")
        if(Bsing<=0) stop("Bsing must be positive.")
        if(Bsing!=round(Bsing)) stop("Bsing must be an integer.")
        if(missing(ncmb)) stop("ncmb is missing.")
        if(ncmb<=0) stop("ncmb must be positive.")
        if(ncmb!=round(ncmb)) stop("ncmb must be an integer.")
        if(ncmb>=dim(D)[1]) stop("ncmb must be smaller than the number of rows in D.")
        if(ncmb<=nsing) stop("ncmb must be higher than nsing.")
        Mseq<-sort(Mseq)
        Mlen<-length(Mseq)
        p<-dim(D)[2]-1
        stabmat<-matrix(numeric((p+1)*Mlen),nrow=Mlen)
        ntrain<-dim(D)[1]
        aggnu<-numeric(p+1)

        for(i in 1:B){
              ind<-sample(1:ntrain,ncmb,replace=FALSE)
              Dcmb<-D[ind,]
              for(m in 1:Mlen){
                  cmbfit<-CMB(D,nsing,Bsing,alpha,singfam,evalfam,sing=sing,M=Mseq[m],m_iter=m_iter,kap=kap,LS=LS,best=best,robagg=robagg,lower=lower)
                  aggnu<-(i-1)/i*aggnu+1/i*cmbfit[[1]]
                  stabmat[m,]<-aggnu
             }
        }
        relev<-rep(list(list()),Mlen)
        for(m in 1:Mlen){
                 relev[[m]]<-which((stabmat[m,]>0)&(stabmat[m,]<1))
        }
        ind<-unique(unlist(relev))
        matplot(Mseq,stabmat[,ind],type='l',yaxt='n',ylim=c(0,1))
        axis(2,at=stabmat[1,ind],las=1,labels=colnames(D)[ind-1])
        axis(4,at=c(0.2,0.4,0.6,0.8),labels=c("0.2","0.4","0.6","0.8"))
        return(list(relev,ind))
}
