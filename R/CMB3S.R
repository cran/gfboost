#' Column Measure Boosting with SingBoost and Stability Selection (CMB-3S)
#'
#' @description{Executes CMB and the loss-based Stability Selection.}
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
#' @param ncmb Number of samples used for \code{CMB}. Integer that must be smaller than the number of samples in \code{Dtrain}.
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
#' @details{See \code{CMB} and \code{CMB.Stabsel}.}
#' @return \item{Final coefficients}{The coefficients corresponding to the optimal stable model as a vector.}
#' \item{Stable column measure}{Aggregated empirical column measure (i.e., selection frequencies) as a vector.}
#' \item{Selected columns}{The column numbers of the variables that form the best stable model as a vector.}
#' \item{Used row measure}{Aggregated empirical row measure (i.e., row weights) as a vector.}
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
#'
#' @examples \donttest{firis<-as.formula(Sepal.Length~.)
#' Xiris<-model.matrix(firis,iris)
#' Diris<-data.frame(Xiris[,-1],iris$Sepal.Length)
#' colnames(Diris)[6]<-"Y"
#' set.seed(19931023)
#' ind<-sample(1:150,120,replace=FALSE)
#' Dtrain<-Diris[ind,]
#' Dvalid<-Diris[-ind,]
#' set.seed(19931023)
#' cmb3s<-CMB3S(Dtrain,nsing=120,Dvalid=Dvalid,ncmb=120,Bsing=1,B=1,alpha=1,singfam=Gaussian()
#' ,evalfam=Gaussian(),sing=FALSE,M=10,m_iter=100,kap=0.1,LS=FALSE,wagg='weights1',
#' gridtype='pigrid',grid=seq(0.8,0.9,1),robagg=FALSE,lower=0,singcoef=TRUE,Mfinal=10)
#' cmb3s$Fin
#' cmb3s$Stab
#' cmb3s$Sel
#' glmres4<-glmboost(Sepal.Length~.,iris[ind,])
#' coef(glmres4)}
#' @examples \donttest{set.seed(19931023)
#' cmb3s1<-CMB3S(Dtrain,nsing=80,Dvalid=Dvalid,ncmb=100,Bsing=10,B=100,alpha=0.5,singfam=Gaussian(),
#' evalfam=Gaussian(),sing=FALSE,M=10,m_iter=100,kap=0.1,LS=FALSE,wagg='weights1',gridtype='pigrid',
#' grid=seq(0.8,0.9,1),robagg=FALSE,lower=0,singcoef=TRUE,Mfinal=10)
#' cmb3s1$Fin
#' cmb3s1$Stab}
#' @examples \donttest{## This will may take around a minute
#' set.seed(19931023)
#' cmb3s2<-CMB3S(Dtrain,nsing=80,Dvalid=Dvalid,ncmb=100,Bsing=10,B=100,alpha=0.5,singfam=Rank(),
#' evalfam=Rank(),sing=TRUE,M=10,m_iter=100,kap=0.1,LS=TRUE,wagg='weights2',gridtype='pigrid',
#' grid=seq(0.8,0.9,1),robagg=FALSE,lower=0,singcoef=TRUE,Mfinal=10)
#' cmb3s2$Fin
#' cmb3s2$Stab}
#' @examples \donttest{set.seed(19931023)
#' cmb3s3<-CMB3S(Dtrain,nsing=80,Dvalid=Dvalid,ncmb=100,Bsing=10,B=100,alpha=0.5,singfam=Huber(),
#' evalfam=Huber(),sing=FALSE,M=10,m_iter=100,kap=0.1,LS=FALSE,wagg='weights2',gridtype='pigrid',
#' grid=seq(0.8,0.9,1),robagg=FALSE,lower=0,singcoef=FALSE,Mfinal=10)
#' cmb3s3$Fin
#' cmb3s3$Stab}
#'
CMB3S<-function(Dtrain,nsing,Bsing=1,B=100,alpha=1,singfam=Gaussian(),evalfam=Gaussian(),sing=FALSE,M=10,m_iter=100,kap=0.1,LS=FALSE,best=1,wagg,gridtype,grid,Dvalid,ncmb,robagg=FALSE,lower=0,singcoef=FALSE,Mfinal=10,...){
        if(missing(Dtrain)) stop("No training data argument.")
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

        cmbstab<-CMB.Stabsel(Dtrain=Dtrain,nsing=nsing,Bsing=Bsing,B=B,alpha=alpha,singfam=singfam,evalfam=evalfam,sing=sing,M=M,m_iter=m_iter,kap=kap,LS=LS,best=best,wagg=wagg,gridtype=gridtype,grid=grid,Dvalid=Dvalid,ncmb=ncmb,robagg=robagg,lower=lower,singcoef=singcoef,Mfinal=Mfinal)
        opt.cols<-cmbstab[[1]]
        opt.coeffs<-cmbstab[[2]]
        names(opt.coeffs)[1]<-"Intercept"
        names(opt.coeffs)[-1]<-colnames(Dtrain)[opt.cols]
        aggnu<-cmbstab[[3]]
        aggzeta<-cmbstab[[4]]

        return(list("Final coefficients"=opt.coeffs,"Stable column measure"=aggnu,"Selected columns"=opt.cols,"Used row measure"=aggzeta))
}
