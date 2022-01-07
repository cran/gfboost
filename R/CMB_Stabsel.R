#' Loss-adapted Stability Selection
#'
#' @description{Workhorse function for the Stability Selection variant where either a grid of thresholds or a grid of
#' cardinalities is given so that the Boosting models are evaluated on a validation set according to all elements of
#' the respective grid. The model which performs best is finally selected as stable model.}
#'
#' @import stats
#'
#' @param Dtrain Data matrix. Has to be an \eqn{n \times (p+1)-}dimensional data frame in the format \eqn{(X,Y)}. The \eqn{X-}part must not
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
#' @param Dvalid Validation data for selecting the optimal element of the grid and with it the best corresponding model.
#' @param ncmb Number of samples used for \code{CMB}. Integer that must be smaller than the number of samples in \code{Dtrain} and higher than \code{nsing}.
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
#' @details{The Stability Selection in the packages \code{stabs} and \code{mboost} requires to fix two of three parameters which are
#' the per-family error rate, the threshold and the number of variables which have to be selected in each model. Our
#' Stability Selection is based on another idea. We also train Boosting models on subsamples but we use a validation
#' step to determine the size of the optimal model. More precisely, if  \code{'pigrid'} is used as \code{gridtype}, the corresponding
#' stable models for each threshold are computed by selecting all variables whose aggregated selection frequency exceeds
#' the threshold. Then, these candidate stable models are validated according to the target loss function (inserted
#' through \code{evalfam}) and the optimal one is finally selected. If \code{'qgrid'} is used as \code{gridtype}, a vector of positive
#' integers has to be entered instead of a vector of thresholds. The candidate stable models then consist of the best
#' variables ordered by their aggregated selection frequencies, respectively. The validation step is the same.}
#' @return \item{colind.opt}{The column numbers of the variables that form the best stable model as a vector.}
#' \item{coeff.opt}{The coefficients corresponding to the optimal stable model as a vector.}
#' \item{aggnu}{Aggregated empirical column measure (i.e., selection frequencies) as a vector.}
#' \item{aggzeta}{Aggregated empirical row measure (i.e., row weights) as a vector.}
#' @export
#'
#' @references{Werner, T., Gradient-Free Gradient Boosting, PhD Thesis, Carl von Ossietzky University Oldenburg, 2020}
#' @references{T. Hothorn, P. Bühlmann, T. Kneib, M. Schmid, and B. Hofner. mboost: Model-Based
#' Boosting, 2017}
#' @references{B. Hofner and T. Hothorn. stabs: Stability Selection with Error Control, 2017.}
#' @references{B. Hofner, L. Boccuto, and M. Göker. Controlling false discoveries in high-dimensional
#' situations: Boosting with stability selection. BMC Bioinformatics, 16(1):144, 2015.}
#' @references{N. Meinshausen and P. Bühlmann. Stability selection. Journal of the Royal Statistical Society:
#' Series B (Statistical Methodology), 72(4):417–473, 2010.}

CMB.Stabsel<-function(Dtrain,nsing,Bsing=1,B=100,alpha=1,singfam=Gaussian(),evalfam=Gaussian(),sing=FALSE,M=10,m_iter=100,kap=0.1,LS=FALSE,best=1,wagg,gridtype,grid,Dvalid,ncmb,robagg=FALSE,lower=0,singcoef=FALSE,Mfinal,...){

      if(missing(Dtrain)) stop("No training data argument.")
      if(missing(nsing)) stop("nsing or ind is missing.")
      if(M!=round(M)) stop("M must be an integer.")
      if(M<=0) stop("M must be positive.")
      if(m_iter!=round(m_iter)) stop("m_iter must be an integer.")
      if(M>m_iter) stop("M must be at most as large as m_iter.")
      if(m_iter<=0) stop("m_iter must be positive.")
      if(kap<=0|kap>1) stop("Learning rate kap must be in ]0,1].")
      if(best<=0|best>1) stop("Parameter best must be in ]0,1].")
      if(nsing<=0) stop("nsing must be positive.")
      if(nsing!=round(nsing)) stop("nsing must be an integer.")
      if(nsing>dim(Dtrain)[1]) stop("nsing must not be higher than the number of rows in Dtrain.")
      if(Bsing<=0) stop("Bsing must be positive.")
      if(Bsing!=round(Bsing)) stop("Bsing must be an integer.")
      if(B<=0) stop("B must be positive.")
      if(B!=round(B)) stop("B must be an integer.")
      if(missing(grid)|missing(gridtype)) stop("grid or gridtype is missing.")
      if(missing(Dvalid)) stop("A validation data set (Dvalid) is missing.")
      if(missing(ncmb)) stop("ncmb is missing.")
      if(ncmb<=0) stop("ncmb must be positive.")
      if(ncmb!=round(ncmb)) stop("ncmb must be an integer.")
      if(ncmb>dim(Dtrain)[1]) stop("ncmb must not be higher than the number of rows in Dtrain.")
      if(ncmb<nsing) stop("ncmb must be at least as high as nsing.")
      if(Mfinal!=round(Mfinal)) stop("Mfinal must be an integer.")
      if(Mfinal<=0) stop("Mfinal must be positive.")

      p<-dim(Dtrain)[2]-1
      ntrain<-dim(Dtrain)[1]
      nvalid<-dim(Dvalid)[1]
      aggnu<-numeric(p+1)
      aggzeta<-numeric(ntrain)

       ## length(CVind.inner)=B
      for(i in 1:B){
            ind<-sample(1:ntrain,ncmb,replace=FALSE)
            Dcmb<-Dtrain[ind,]
            cmb<-CMB(D=Dcmb,nsing=nsing,Bsing=Bsing,alpha=alpha,singfam=singfam,evalfam=evalfam,sing=sing,M=M,m_iter=m_iter,kap=kap,LS=LS,best=best,wagg=wagg,robagg=robagg,lower=lower)
            aggnu<-(i-1)/i*aggnu+1/i*cmb[[1]]
            aggzeta<-(i-1)/i*aggzeta
            aggzeta[ind]<-aggzeta[ind]+1/i*cmb[[4]]
      }

     len<-length(grid)

     loss<-numeric(len)

     for(k in 1:len){
            if(gridtype=='qgrid'){
                  ## q columns with the highest selection frequencies
                  colind<-which(aggnu>=sort(aggnu,decreasing=TRUE)[grid[k]])-1
            }
            if(gridtype=='pigrid'){
                  ## all columns with selection frequencies larger than the threshold
                  colind<-which(aggnu>=grid[k])-1
            }
           ## catch the case that no column is selected. An empty model is very problematic so that we set the loss
           ## to infinity to guarantee that such a model will never be the best model (unless all models are empty)
           if(length(colind)==0){
                 loss[k]<-Inf
           }
           else if(min(colind)==0){
                 loss[k]<-Inf
           }
           else{
                 if(singcoef==TRUE){
                      res<-singboost(Dtrain[,c(colind,p+1)],M=M,m_iter=m_iter,kap=kap,best=best,LS=LS,singfamily=singfam)
                      pred<-as.matrix(cbind(rep(1,nvalid),Dvalid[,colind]))%*%res$Coeff
                      if(evalfam@name!='Localized ranking Family'){
                             loss[k]<-evalfam@risk(Dvalid$Y,pred)
                      }
                      if(evalfam@name=='Localized ranking Family'){
                             loss[k]<-LocRank(ceiling(best*nvalid))@risk(Dvalid$Y,pred)
                      }
                  }
                 if(singcoef==FALSE){
                       reslm<-lm(Y~.,Dtrain[,c(colind,p+1)],weights=aggzeta)
                       if(evalfam@name!='Localized ranking Family'){
                              loss[k]<-evalfam@risk(Dvalid$Y,predict(reslm,Dvalid))
                       }
                       if(evalfam@name=='Localized ranking Family'){
                              loss[k]<-LocRank(ceiling(best*nvalid))@risk(Dvalid$Y,predict(reslm,Dvalid))
                       }
                 }
           }
      }

      if(min(loss)==Inf){
            warning("Infinite loss")
            colind.opt<-which.max(aggnu[-1])
            if(singcoef==FALSE){
                 coeff.opt<-coef(lm(Y~.,Dtrain[,c(colind.opt,p+1)],weights=aggzeta))
            }
            if(singcoef==TRUE){
                 coeff.opt<-singboost(Dtrain[,c(colind.opt,p+1)],M=Mfinal,m_iter=m_iter,kap=kap,best=best,LS=LS,singfamily=singfam)$Coeff
            }
      }

      if(min(loss)!=Inf){
           k.opt<-which.min(loss)
           ## find again the optimal columns according to k.opt
           if(gridtype=='qgrid'){
                colind.opt<-which(aggnu>=sort(aggnu,decreasing=TRUE)[grid[k.opt]])-1
           }
           if(gridtype=='pigrid'){
                colind.opt<-which(aggnu>=grid[k.opt])-1
           }
           ## weighted least squares model where the rows are weighted according to the aggregated row measure aggzeta
           if(singcoef==FALSE){
                coeff.opt<-coef(lm(Y~.,Dtrain[,c(colind.opt,p+1)],weights=aggzeta))
           }
           if(singcoef==TRUE){
                coeff.opt<-singboost(Dtrain[,c(colind.opt,p+1)],M=Mfinal,m_iter=m_iter,kap=kap,best=best,LS=LS,singfamily=singfam)$Coeff
           }
      }


    return(list(colind.opt,coeff.opt,aggnu,aggzeta))
}
