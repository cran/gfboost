#' Plot function for the SingBoost coefficient paths
#'
#' @import stats
#' @import graphics
#'
#' @param mod singboost object.
#' @param M An integer between 2 and \code{m_iter}. Indicates that in every M-th iteration, a singular iteration will be
#' performed. Default is 10.
#' @param m_iter Number of SingBoost iterations. Default is 100.
#' @param subnames Use it only if the variable names are of the form ''letter plus number''. Better just ignore it.
#'
#' @return Nothing. Plots SingBoost coefficient paths
#' @export
#'
#' @examples {glmres<-glmboost(Sepal.Length~.,iris)
#' glmres
#' attributes(varimp(glmres))$self
#' attributes(varimp(glmres))$var
#' firis<-as.formula(Sepal.Length~.)
#' Xiris<-model.matrix(firis,iris)
#' Diris<-data.frame(Xiris[,-1],iris$Sepal.Length)
#' plot(glmres)
#' singpath<-path.singboost(Diris)
#' singboost.plot(singpath,10,100,subnames=FALSE)}
#'
singboost.plot<-function(mod,M,m_iter,subnames=FALSE){
    cp<-mod[[5]]
    unq<-unique(cp[2,])
    pathlen<-length(unique(cp[2,]))
    pathlist<-rep(list(list()),pathlen)
    if(subnames==TRUE){
       pathnames<-unq[order(as.numeric(substring(unq,first=2)))]
    }
    else{pathnames<-unq}
    names(pathlist)<-pathnames
    for(i in 1:pathlen){
             Z<-numeric(m_iter)
             Z[which(cp[2,]==pathnames[i])]<-as.numeric(cp[1,][which(cp[2,]==pathnames[i])])
             pathlist[[i]]<-Z
             }
    pathplot<-lapply(pathlist,cumsum)
    pathmat<-matrix(unlist(pathplot),ncol=pathlen)
    matplot(pathmat,type='l',main=paste("Singboost coefficient paths \n",  mod[[6]]@name),
    ylab="Coefficients", xlab="Number of iterations")
    #points(cumsum(res[[4]]),type='l')
    axis(4,at=pathmat[m_iter,],labels=pathnames,las=1)
    #axis(4,at=res[[4]][m_iter],labels="Intercept",las=1)
    abline(h = 0, lty = 1, col = "lightgray")
    abline(v=seq(1,m_iter-M+1,by=M),col='lightgray')
}
