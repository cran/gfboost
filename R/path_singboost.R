#' Coefficient paths for SingBoost
#'
#' @description{Runs SingBoost but saves the coefficients paths. If no coefficient path plot is needed, just use
#' \code{singboost}.}
#'
#' @import mboost
#' @param D Data matrix. Has to be an \eqn{n \times (p+1)-}dimensional data frame in the format \eqn{(X,Y)}. The \eqn{X-}part must not
#' contain an intercept column containing only ones since this column will be added automatically.
#' @param M An integer between 2 and \code{m_iter}. Indicates that in every \eqn{M-}th iteration, a singular iteration will be
#' performed. Default is 10.
#' @param m_iter Number of SingBoost iterations. Default is 100.
#' @param kap Learning rate (step size). Must be a real number in \eqn{]0,1]}. Default is 0.1 It is recommended to use
#' a value smaller than 0.5.
#' @param singfamily A Boosting family corresponding to the target loss function. See .\code{mboost} for families
#' corresponding to standard loss functions. May also use the loss functions for ranking losses provided in this
#' package. Default is \code{Gaussian()} for which SingBoost is just standard \eqn{L_2-}Boosting.
#' @param best Needed in the case of localized ranking. The parameter \code{K} of the localized ranking loss will be
#' computed by \eqn{best \cdot n} (rounded to the next larger integer). Warning: If a parameter \code{K} is inserted into the
#' \code{LocRank} family, it will be ignored when executing SingBoost.
#' @param LS If a \code{singfamily} object that is already provided by \code{mboost} is used, the respective Boosting algorithm
#' will be performed in the singular iterations if \code{Ls} is set to \code{TRUE}. Default is \code{FALSE}.
#'
#'
#' @return \item{Selected variables}{Names of the selected variables.}
#' \item{Coefficients}{The selected coefficients as an \eqn{(p+1)-}dimensional vector (i.e., including the zeroes).}
#' \item{Freqs}{Selection frequencies and a matrix for intercept and coefficient paths, respectively.}
#' \item{Intercept path}{The intercept path as an \eqn{m_{iter}-}dimensional vector.}
#' \item{Coefficient path}{The coefficient paths as a \eqn{2 \cdot m_{iter} \times 2-}dimensional matrix.}
#' @export
#'
#'
 path.singboost<-function(D,M=10,m_iter=100,kap=0.1,singfamily=Gaussian(),best=1,LS=FALSE){
     if(M!=round(M)) stop("M must be an integer.")
     if(M<=0) stop("M must be positive.")
     if(M>m_iter) stop("M must be at most as large as m_iter.")
     if(m_iter!=round(m_iter)) stop("m_iter must be an integer.")
     if(m_iter<=0) stop("m_iter must be positive.")
     if(kap<=0|kap>1) stop("Learning rate kap must be in ]0,1].")
     if(best<=0|best>1) stop("Parameter best must be in ]0,1].")
     p<-dim(D)[2]-1
     ## vector for coefficients
     coeff<-numeric(p+1)
     ## matrix for coefficient paths
     coeffpath<-matrix(numeric(2*m_iter),nrow=2)
     ## vector for intercept path
     interpath<-numeric(m_iter)
     runs<-floor(m_iter/M)
     ## vector for selection frequencies of the variables (including the intercept)
     occ<-numeric(p+1)
     n<-dim(D)[1]
     x<-D[,1:p]
     X<-cbind(rep(1,n),x)
     if (p==1) colnames(X)[2]<-colnames(D)[1]
     colnames(D)[p+1]<-"Y"
     y<-D$Y
     r<-y
     for (i in 1:runs){
           ###################### Singular iteration ######################
           if(singfamily@name%in%c("Ranking family","Localized ranking family")==FALSE){
               if(LS==FALSE){
                  res<-glmboost((Y-mean(Y))~.,data=D,family=singfamily,control=
                           boost_control(mstop=1,nu=kap))
                  ## centering: otherwise, the first (and only in this case)
                  ## iteration would be the computation of the offset
                  coeff[1]<-coeff[1]+mean(r)+coef(res)[1]
                  ## correct the centering from above
                if(length(names(coef(res)))==2){
                     j<-which(colnames(X)==names(coef(res)[2]))
                }
                if(length(names(coef(res)))==1){
                     j<-1
                 }
                 ## For the special case that p=1, x has no colname,
                 ## while the usage of X does not cause such problems
                 if(j>1) coeff[j]<-coeff[j]+coef(res)[2]
                 ## catches the case that j=1, so that the intercept would
                 ## be updated two times
                 occ[j]<-occ[j]+1
                 interpath[M*(i-1)+1]<-mean(r)
                 coeffpath[1,M*(i-1)+1]<-coef(res)[2]
                 coeffpath[2,M*(i-1)+1]<-names(coef(res)[2])
                 r<-r-mean(r)-coef(res)[1]-coef(res)[2]*X[,j]
              }
              if(LS==TRUE){
                  losses<-numeric(p)
                  for(k in 1:p){
                      losses[k]<-singfamily@risk(predict(lm(r~D[,k])),r)
                  }
                  ## for loss functions that satisfy L(y,xbeta)=L(y-xbeta)
                  ## it is sufficient to select the baselearner that models
                  ## the residuals best (''inner loss'')
                  j<-as.numeric(sample(as.character(which.min(losses)),1))
                  res<-lm(r~D[,j])
                  ## Refit the best model to avoid having saved all models

                  coeff[1]<-coeff[1]+kap*coef(res)[1]
                  coeff[j+1]<-coeff[j+1]+kap*coef(res)[2]
                  ## the learning rate must be respected here
                  ## for LS=FALSE and a standard loss function, glmboost did this automatically
                  occ[j+1]<-occ[j+1]+1
                  interpath[M*(i-1)+1]<-mean(r)
                  coeffpath[1,M*(i-1)+1]<-coef(res)[2]
                  coeffpath[2,M*(i-1)+1]<-names(coef(res))[2]
                  r<-r-kap*coef(res)[1]-kap*coef(res)[2]*X[,j+1]
             }
        }
        if(singfamily@name=="Ranking family"){
              losses<-numeric(p)
              for(k in 1:p){
                    losses[k]<-Rank()@risk(y-r+kap*predict(lm(r~D[,k])),y)
              }
             ## For ranking problems, the ''inner loss'' is not sufficient.
             ## Therefore, we have to check the performance of the updated
             ## model (''outer loss'') for each baselearner

             ## the ranking loss is discrete, so ties are possible:
             ## break them randomly by uniform sampling
             j<-as.numeric(sample(as.character(which.min(losses)),1))
             res<-lm(r~D[,j])
             coeff[1]<-coeff[1]+kap*coef(res)[1]
             coeff[j+1]<-coeff[j+1]+kap*coef(res)[2]
             occ[j+1]<-occ[j+1]+1
             interpath[M*(i-1)+1]<-mean(r)
             coeffpath[1,M*(i-1)+1]<-kap*coef(res)[2]
             coeffpath[2,M*(i-1)+1]<-colnames(D)[j]
             r<-r-kap*coef(res)[1]-kap*coef(res)[2]*X[,j+1]
        }
        if(singfamily@name=="Localized ranking family"){
                losses<-numeric(p)
                for(k in 1:p){
                    losses[k]<-LocRank(ceiling(best*n))@risk(y-r+kap*predict(lm(r~D[,k])),y)
                }
                ## the ranking loss is discrete, so ties are possible:
                ## break them randomly by uniform sampling
                j<-as.numeric(sample(as.character(which.min(losses)),1))
                res<-lm(r~D[,j])
                coeff[1]<-coeff[1]+kap*coef(res)[1]
                coeff[j+1]<-coeff[j+1]+kap*coef(res)[2]
                occ[j+1]<-occ[j+1]+1
                interpath[M*(i-1)+1]<-mean(r)
                coeffpath[1,M*(i-1)+1]<-kap*coef(res)[2]
                coeffpath[2,M*(i-1)+1]<-colnames(D)[j]
                r<-r-kap*coef(res)[1]-kap*coef(res)[2]*X[,j+1]
       }
       D$Y<-r
       ## End of singular iteration, get current residuals

       ####################### L2-Iterations ########################
       res<-glmboost((Y-mean(Y))~.,data=D,family=Gaussian(),control=
                    boost_control(mstop=M-1,nu=kap))
       coeff[1]<-coeff[1]+mean(r)+coef(res)[1]
       j<-which(colnames(X)%in%c(names(coef(res))[-1]))
        coeff[j]<-coeff[j]+coef(res)[-1]
        ## as before: the (M-1) Boosting iterations and the
        ## computed coefficients already respect the learning rate
        occ<-occ+attributes(varimp(res))$selfreqs*(M-1)
        ## get the selection frequencies with varimp
        u<-coef(res,aggregate='none')
        N<-matrix(numeric(2*(M-1)),nrow=2)
        for(o in 2:length(u)){
              jnd<-which(u[[o]]!=0)
              N[1,jnd]<-u[[o]][jnd]
              N[2,jnd]<-names(u)[o]
        }
        interpath[(M*(i-1)+2):(M*i)]<-c(u[[1]])
        coeffpath[1,(M*(i-1)+2):(M*i)]<-N[1,]
        coeffpath[2,(M*(i-1)+2):(M*i)]<-N[2,]
        r<-y-as.matrix(X[,occ!=0])%*%coeff[occ!=0]-coeff[1]
        D$Y<-r
        ## End of standard iterations, get current residuals
     }
     ## get the positions of the relevant variables detected by Singboost
     selind<-which(occ!=0)


     ## get the names of the relevant variables and order them by their selection frequencies
     selvars<-paste(c("Intercept",colnames(D))[which(occ!=0)][order(occ[occ!=0],decreasing=TRUE)],collapse='>=')

     ## Sel: Names of the selected variables
     ## Coef: The selected coefficients as an (p+1)-dimensional vector (i.e., including the zeroes)
     ## Freqs: The relative selection frequencies for each column (i.e., the empirical column measure

     return(list("Selected variables"=selvars,"Coefficients"=coeff,"Freqs"=occ/(ceiling(m_iter/M)*M),"Intercept path"=interpath,
  "Coefficient path"=coeffpath,singfamily))
}
