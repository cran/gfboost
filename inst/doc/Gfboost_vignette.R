## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(gfboost)

## ----eval=TRUE----------------------------------------------------------------
y<-c(-3, 10.3,-8, 12, 14,-0.5, 29,-1.1,-5.7, 119)
yhat<-c(0.02, 0.6, 0.1, 0.47, 0.82, 0.04, 0.77, 0.09, 0.01, 0.79)
Rank()@risk(y,yhat)

## ----eval=TRUE----------------------------------------------------------------
glmres<-glmboost(Sepal.Length~.,iris)
glmres
attributes(varimp(glmres))$self
attributes(varimp(glmres))$var
firis<-as.formula(Sepal.Length~.)
Xiris<-model.matrix(firis,iris)
Diris<-data.frame(Xiris[,-1],iris$Sepal.Length)
colnames(Diris)[6]<-"Y"
coef(glmboost(Xiris,iris$Sepal.Length))
singboost(Diris)
singboost(Diris,LS=TRUE)

## ----eval=TRUE----------------------------------------------------------------
glmres2<-glmboost(Sepal.Length~Petal.Length+Sepal.Width:Species,iris)
finter<-as.formula(Sepal.Length~Petal.Length+Sepal.Width:Species-1)
Xinter<-model.matrix(finter,iris)
Dinter<-data.frame(Xinter,iris$Sepal.Length)
singboost(Dinter)
coef(glmres2)

## ----eval=TRUE----------------------------------------------------------------
glmquant<-glmboost(Sepal.Length~.,iris,family=QuantReg(tau=0.75))
coef(glmquant)
attributes(varimp(glmquant))$self
singboost(Diris,singfamily=QuantReg(tau=0.75),LS=TRUE)
singboost(Diris,singfamily=QuantReg(tau=0.75),LS=TRUE,M=2)
singboost(Diris,singfamily=Rank(),LS=TRUE)

## ----eval=TRUE----------------------------------------------------------------
singpath<-path.singboost(Diris)
singboost.plot(singpath,10,100)

## ----eval=TRUE----------------------------------------------------------------
set.seed(19931023)
cmb1<-CMB(Diris,nsing=100,Bsing=50,alpha=0.8,singfam=Rank(),
evalfam=Rank(),sing=TRUE,M=10,m_iter=100,
kap=0.1,LS=TRUE,wagg='weights1',robagg=FALSE,lower=0)
cmb1

## ----eval=TRUE----------------------------------------------------------------
set.seed(19931023)
ind<-sample(1:150,120,replace=FALSE)
Dtrain<-Diris[ind,]
Dvalid<-Diris[-ind,]
set.seed(19931023)
cmb3s1<-CMB3S(Dtrain,nsing=80,Dvalid=Dvalid,ncmb=80,Bsing=1,B=100,alpha=0.5,singfam=Gaussian(),
evalfam=Gaussian(),sing=FALSE,M=10,m_iter=100,kap=0.1,LS=FALSE,wagg='weights1',gridtype='pigrid',
grid=seq(0.8,0.9,1),robagg=FALSE,lower=0,singcoef=TRUE,Mfinal=10)
cmb3s1$Fin
cmb3s1$Stab

## ----eval=TRUE----------------------------------------------------------------
set.seed(19931023)
ind<-sample(1:150,120,replace=FALSE)
Dtrain<-Diris[ind,]
Dvalid<-Diris[-ind,]
set.seed(19931023)
cmb3s2<-CMB3S(Dtrain,nsing=80,Dvalid=Dvalid,ncmb=80,Bsing=1,B=100,alpha=0.5,singfam=Gaussian(),
evalfam=Gaussian(),sing=FALSE,M=10,m_iter=100,kap=0.1,LS=FALSE,wagg='weights1',gridtype='qgrid',
grid=seq(1,2,3),robagg=FALSE,lower=0,singcoef=TRUE,Mfinal=10)
cmb3s2$Fin
cmb3s2$Stab

