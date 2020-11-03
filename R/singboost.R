#' SingBoost Boosting method
#'
#' @import mboost
#' @description{SingBoost is a Boosting method that can deal with complicated loss functions that do not allow for
#' a gradient. SingBoost is based on L2-Boosting in its current implementation.}
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
#' @details{Gradient Boosting algorithms require convexity and differentiability of the underlying loss function.
#' SingBoost is a Boosting algorithm based on \eqn{L_2-}Boosting that allows for complicated loss functions that do not
#' need to satisfy these requirements. In fact, SingBoost alternates between standard \eqn{L_2-}Boosting iterations and
#' singular iterations where essentially an empirical gradient step is executed in the sense that the baselearner
#' that performs best, evaluated in the complicated loss, is selected in the respective iteration. The implementation
#' is based on \code{glmboost} from the package \code{mboost} and using the \eqn{L_2-}loss in the singular iterations returns exactly the
#' same coefficients as \eqn{L_2-}Boosting.}
#'
#' @return \item{Selected variables}{Names of the selected variables.}
#' \item{Coefficients}{The selected coefficients as an \eqn{(p+1)-}dimensional vector (i.e., including the zeroes).}
#' \item{Freqs}{Selection frequencies and a matrix for intercept and coefficient paths, respectively.}
#' \item{VarCoef}{Vector of the non-zero coefficients.}
#' @export
#'
#' @references{Werner, T., Gradient-Free Gradient Boosting, PhD Thesis, Carl von Ossietzky University Oldenburg, 2020}
#' @references{P. Bühlmann and B. Yu. Boosting with the l2 loss: Regression and Classification. Journal
#' of the American Statistical Association, 98(462):324–339, 2003}
#' @references{T. Hothorn, P. Bühlmann, T. Kneib, M. Schmid, and B. Hofner. mboost: Model-Based
#' Boosting, 2017}
#'
#' @examples {glmres<-glmboost(Sepal.Length~.,iris)
#' glmres
#' attributes(varimp(glmres))$self
#' attributes(varimp(glmres))$var
#' firis<-as.formula(Sepal.Length~.)
#' Xiris<-model.matrix(firis,iris)
#' Diris<-data.frame(Xiris[,-1],iris$Sepal.Length)
#' colnames(Diris)[6]<-"Y"
#' coef(glmboost(Xiris,iris$Sepal.Length))
#' singboost(Diris)
#' singboost(Diris,LS=TRUE)}
#' @examples {glmres2<-glmboost(Sepal.Length~Petal.Length+Sepal.Width:Species,iris)
#' finter<-as.formula(Sepal.Length~Petal.Length+Sepal.Width:Species-1)
#' Xinter<-model.matrix(finter,iris)
#' Dinter<-data.frame(Xinter,iris$Sepal.Length)
#' singboost(Dinter)
#' coef(glmres2)}
#' @examples {glmres3<-glmboost(Xiris,iris$Sepal.Length,control=boost_control(mstop=250,nu=0.05))
#' coef(glmres3)
#' attributes(varimp(glmres3))$self
#' singboost(Diris,m_iter=250,kap=0.05)
#' singboost(Diris,LS=TRUE,m_iter=250,kap=0.05)}
#' @examples {glmquant<-glmboost(Sepal.Length~.,iris,family=QuantReg(tau=0.75))
#' coef(glmquant)
#' attributes(varimp(glmquant))$self
#' singboost(Diris,singfamily=QuantReg(tau=0.75),LS=TRUE)
#' singboost(Diris,singfamily=QuantReg(tau=0.75),LS=TRUE,M=2)}
#' @examples {singboost(Diris,singfamily=Rank(),LS=TRUE)
#' singboost(Diris,singfamily=Rank(),LS=TRUE,M=2)}


singboost<-function(D,M=10,m_iter=100,kap=0.1,singfamily=Gaussian(),best=1,LS=FALSE){

    if(missing(D)) stop("No data argument.")
    if(M!=round(M)) stop("M must be an integer.")
    if(M<=0) stop("M must be positive.")
    if(m_iter!=round(m_iter)) stop("m_iter must be an integer.")
    if(M>m_iter) stop("M must be at most as large as m_iter.")
    if(m_iter<=0) stop("m_iter must be positive.")
    if(kap<=0|kap>1) stop("Learning rate kap must be in ]0,1].")
    if(best<=0|best>1) stop("Parameter best must be in ]0,1].")
    p<-dim(D)[2]-1
    ## vector for coefficients
    coeff<-numeric(p+1)
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
    r<-y-as.matrix(X[,occ!=0])%*%coeff[occ!=0]-coeff[1]
    D$Y<-r
    ## End of standard iterations, get current residuals
   }
   ## get the positions of the relevant variables detected by Singboost
   selind<-which(occ!=0)

   ## get the names of the relevant variables and order them by their selection frequencies
   selvars<-paste(c("Intercept",colnames(D))[which(occ!=0)][order(occ[occ!=0],decreasing=T)],collapse='>=')
   varcoef<-coeff[coeff!=0]
   names(varcoef)<-c("",colnames(D)[which(occ!=0)-1])
   names(varcoef)[1]<-"Intercept"

   ## Sel: Names of the selected variables
   ## Coef: The selected coefficients as an (p+1)-dimensional vector (i.e., including the zeroes)
   ## Freqs: The relative selection frequencies for each column (i.e., the empirical column measure

  return(list("Selected variables"=selvars,"Coefficients"=coeff,"Freqs"=occ/(ceiling(m_iter/M)*M),"VarCoef"=varcoef))
}
