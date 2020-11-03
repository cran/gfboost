#' CMB aggregation function
#'
#' @description{Aggregates the selection frequencies of multiple SingBoost models. May be used with caution since
#' there are not yet recommendations about good hyperparameters.}
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
#' @param robagg Optional. If setting \code{robagg=TRUE}, the best SingBoost models are ignored when executing the
#' aggregation to avoid inlier effects. Only reasonable in combination with \code{lower}.
#' @param lower Optional argument. Only reasonable when setting \code{robagg=TRUE}. \code{lower} is a real number in \eqn{[0,1[} (a rather
#' small number is recommended) and indicates that the aggregation ignores the SingBoost models with the best
#' performances to avoid possible inlier effects.
#' @param ... Optional further arguments
#'
#' @details{SingBoost is designed to detect variables that standard Boosting procedures may not but which may be
#' relevant w.r.t. the target loss function. However, one may try to stabilize this ''singular part'' of the
#' column measure by aggregating several SingBoost models in the sense that they are evaluated on a validation set
#' and that the selection frequencies are averaged, maybe in a weighted manner according to the validation losses.
#' Warning: This procedure does not replace a Stability Selection!}
#'
#' @return \item{Column measure}{Aggregated column measure as \eqn{(p+1)-}dimensional vector.}
#'  \item{Selected variables}{Names of the variables with positive aggregated column measure.}
#'  \item{Variables names}{Names of all variables including the intercept.}
#'  \item{Row measure}{Aggregated row measure as \eqn{n-}dimensional vector.}
#' @export
#'
#' @references{Werner, T., Gradient-Free Gradient Boosting, PhD Thesis, Carl von Ossietzky University Oldenburg, 2020}
#'
#' @examples \donttest{firis<-as.formula(Sepal.Length~.)
#' Xiris<-model.matrix(firis,iris)
#' Diris<-data.frame(Xiris[,-1],iris$Sepal.Length)
#' colnames(Diris)[6]<-"Y"
#' set.seed(19931023)
#' cmb1<-CMB(Diris,nsing=100,Bsing=50,alpha=0.8,singfam=Rank(),
#' evalfam=Rank(),sing=TRUE,M=10,m_iter=100,
#' kap=0.1,LS=TRUE,wagg='weights1',robagg=FALSE,lower=0)
#' cmb1
#' set.seed(19931023)
#' cmb2<-CMB(Diris,nsing=100,Bsing=50,alpha=0.8,singfam=Rank(),
#' evalfam=Rank(),sing=TRUE,M=2,m_iter=100,
#' kap=0.1,LS=TRUE,wagg='weights1',robagg=FALSE,lower=0)
#' cmb2[[1]]
#' set.seed(19931023)
#' cmb3<-CMB(Diris,nsing=100,Bsing=50,alpha=0.8,singfam=Rank(),
#' evalfam=Rank(),sing=TRUE,M=10,m_iter=100,
#' kap=0.1,LS=TRUE,wagg='weights2',robagg=FALSE,lower=0)
#' cmb3[[1]]}

CMB<-function(D,nsing,Bsing=1,alpha=1,singfam=Gaussian(),evalfam=Gaussian(),sing=FALSE,M=10,m_iter=100,kap=0.1,LS=FALSE,best=1,wagg,robagg=FALSE,lower=0,...){

     if(missing(D)) stop("No data argument.")
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
     if(missing(wagg)) wagg<-"weights1"
     p<-dim(D)[2]-1
     ncmb<-dim(D)[1]

     ind<-sapply(1:Bsing,function(i) sample(1:ncmb,nsing,replace=FALSE))

     ## Step (LS), (Sing), (CoM)

     res.L2<-RejStep(D,nsing=nsing,Bsing=Bsing,ind=ind,sing=sing,singfam=singfam,evalfam=evalfam,M=M,m_iter=m_iter,kap=kap,LS=LS,best=best)
     loss<-res.L2$loss
     occ<-res.L2$occ
     ## Choose the best ceil(alpha*Bsing) Singboost models according to the out-of-sample loss (here, only the number
     ## of the model is saved)
     ind.best<-which(loss<=sort(loss)[ceiling(alpha*Bsing)])
     ## Step (W-AGG)

     occagg<-numeric(p+1)
     zetaagg<-numeric(ncmb)
     if(wagg=="weights1"){
         for(i in ind.best){
                occagg<-occagg+1*(occ[[i]]>0)
                zetaagg[ind[,i]]<-zetaagg[ind[,i]]+1
         }
        ## standardization to [0,1]
        occagg<-occagg/length(ind.best)
        zetaagg<-zetaagg/length(ind.best)
      }
     if(wagg=="weights2"){
         if(robagg==TRUE){
             loss[1:floor(Bsing*lower)]<-loss[1:ceiling(Bsing*lower)]
         }
         wsum<-sum(loss[ind.best])
         for(i in ind.best){
                w<-loss[i]/wsum
                occagg<-occagg+w*1*(occ[[i]]>0)
                zetaagg[ind[,i]]<-zetaagg[ind[,i]]+rep(w,nsing)
                }
          ## no standardization necessary since the weights automatically sum up to 1
     }
     ## aggregated column measure represented by the vector occagg of selection frequencies
     ## aggregated row measure represented by the vector zetaagg of weights


     finalvars<-paste(c("Intercept",colnames(D))[occagg!=0][order(occagg[occagg!=0],decreasing=TRUE)],collapse='>=')
     ## Ordered variable names according to the variables with positive column measure

     ## Output
     ## occagg: aggregated column measure as (p+1)-dimensional vector
     ## finalvars: Names of the variables with positive aggregated column measure
     ## $3: Names of all variables including the intercept
     ## zetaagg: Aggregated row measure as n-dimensional vector
     return(list("Column measure"=occagg,"Selected variables"=finalvars,"Variable names"=c("Intercept",colnames(D)),"Row measure"=zetaagg))
}
