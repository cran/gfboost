#' CMB validation step
#'
#' @description{Validation step to combine different SingBoost models.}
#'
#' @import stats
#'
#' @param D Data matrix. Has to be an \eqn{n \times (p+1)-}dimensional data frame in the format \eqn{(X,Y)}. The \eqn{X-}part must not
#' contain an intercept column containing only ones since this column will be added automatically.
#' @param nsing Number of observations (rows) used for the SingBoost submodels.
#' @param Bsing Number of subsamples based on which the SingBoost models are validated. Default is 1. Not to confuse with parameter \code{B} for the Stability Selection.
#' @param ind Vector with indices for dividing the data set into training and validation data.
#' @param sing If \code{sing=FALSE} and the \code{singfam} family is a standard Boosting family that is contained in the package
#' \code{mboost}, the CMB aggregation procedure is executed for the corresponding standard Boosting models.
#' @param singfam A SingBoost family. The SingBoost models are trained based on the corresponding loss function. Default is \code{Gaussian()} (squared loss).
#' @param evalfam A SingBoost family. The SingBoost models are validated according to the corresponding loss function. Default is \code{Gaussian()} (squared loss).
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
#'
#' @details{Divides the data set into a training and a validation set. The SingBoost models are computed on the training
#' set and evaluated on the validation set based on the loss function corresponding to the selected Boosting family.}
#'
#' @return \item{loss}{Vector of validation losses.}
#' \item{occ}{Selection frequencies for each Boosting model.}
#'
#' @references{Werner, T., Gradient-Free Gradient Boosting, PhD Thesis, Carl von Ossietzky University Oldenburg, 2020}
#'
RejStep<-function(D,nsing,Bsing=1,ind,sing=FALSE,singfam=Gaussian(),evalfam=Gaussian(),M=10,m_iter=100,kap=0.1,LS=FALSE,best=1){

      if(missing(D)) stop("No data argument.")
      if(missing(nsing)|missing(ind)) stop("nsing or ind is missing.")
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
      ncmb<-dim(D)[1]
      p<-dim(D)[2]-1
      loss<-numeric(Bsing)
      occ<-rep(list(list()),Bsing)

      for(k in 1:Bsing){
           Dsing<-D[ind[,k],]
           Dsingtest<-D[-ind[,k],]

          ## K will be needed two times: For the training data in SingBoost and for the test loss when finding the best
          ## Singboost models. Both requires a different K. The first one will be chosen automatically in singboost, the
          ## second one will be, disregarding the input, computed user-friendly according to the parameter best and the
          ## number of observations in the validation data


          ## Step (GLMBOOST)
          ## CMB for standard Boosting families
          if(sing==FALSE){
          if(singfam@name=='Localized ranking family'){
              res<-glmboost(Y~.,Dsing,family=LocRank(ceiling(best*nsing)),control=boost_control(nu=kap,mstop=m_iter))
          }
          if(singfam@name!='Localized ranking family'){
              res<-glmboost(Y~.,Dsing,family=singfam,control=boost_control(nu=kap,mstop=m_iter))
          }
          yhat<-predict(res,Dsingtest)
          occ[[k]]<-attributes(varimp(res))$selfreqs
          }
          ## Step (SING)
          ## CMB for non-standard Boosting families where singular iterations are executed
          if(sing==TRUE){
              res<-singboost(Dsing,M=M,m_iter=m_iter,kap=kap,singfamily=singfam,best=best,LS=LS)
              yhat<-as.matrix(cbind(rep(1,ncmb-nsing),Dsingtest[,1:p]))%*%as.vector(res$Coef)
              occ[[k]]<-res$Freqs
          }
          ## Out-of-sample loss
          if(evalfam@name=='Localized ranking family'){
             loss[k]<-LocRank(ceiling(best*(ncmb-nsing)))@risk(Dsingtest$Y,yhat)
          }
         if(evalfam@name!='Localized ranking family'){
             loss[k]<-evalfam@risk(Dsingtest$Y,yhat)
          }
     }

     ## Step (CoM)
     ## loss: Vector of validation losses
     ## occ: List of selection frequencies for each Singboost model
     return(list("loss"=loss,"occ"=occ))
}
