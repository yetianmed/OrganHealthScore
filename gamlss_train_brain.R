rm(list = ls())
# set up your working directory using setwd("xxx")
# install the following R packages if you don't have using install.packages (xxx)
# you may not need all of them...
library(ggpointdensity)
library(cowplot)
library(patchwork)
library(qqplotr)
library(R.matlab)
library(Hmisc)
library(rstatix)
library(ggpubr)
library(scoring)
library(gamlss)
library(caret)

# Read fake data
M<-readMat("fake_data.mat")

# Run 10-fold CV to compute z-scores for HC individuals 
# 10-fold CV
nfolds<-10
zscores<-matrix(data=NA,nrow=nrow(M$data),ncol=1)
conv<-matrix(data=0,nrow=1,ncol=nfolds)
bic_val<-matrix(data=0,nrow=1,ncol=nfolds)
data<-data.frame(cbind(M$age,M$sex,M$site,M$data))
names(data)<-c("age","sex","site","phenotype")

# loop over 10 folds
folds_idx<-createFolds(c(1:nrow(data)),k=nfolds)
scores<-matrix(data=NA,nrow=nrow(data),1)
  
for (k in 1:nfolds){
  
  idx_test<-folds_idx[[k]]
  idx_train<-setdiff(c(1:nrow(data)),idx_test)
  
  data_test<-data[idx_test,]
  data_train<-data[idx_train,]
  
  # Model fitting
  con<-gamlss.control(n.cyc=100)
  mdl<-gamlss(phenotype~fp(age,npoly=2) + sex + random(as.factor(site)),
              sigma.fo=~fp(age,npoly=2),
              nu.fo=~1,
              nu.tau=~1,
              family=BCT,
              data=data_train,
              control=con)
  if (mdl$converged)
    print("Model converged")
  else
    print("Model not converged")
  
  conv[k]<-mdl$converged
  bic_val[k]<-mdl$sbc
  
  # Predict new data
  # Compute z-score #
  # https://rdrr.io/cran/gamlss.dist/man/BCt.html
  # BCT
  
  params<-predictAll(mdl,data=data_train,
                     newdata=data_test,output='matrix',type="response",
                     y.value="median",what=c("mu", "sigma", "nu","tau"))
  quantiles <- pBCT(params[,1],mu=params[,2],sigma=params[,3],nu=params[,4],tau=params[,5])
  
  # randomized residuals are quantiles transformed to corresponding
  # z-scores of normal distribution
  z_randomized <- qnorm(quantiles, mean = 0, sd = 1)
  qplot(z_randomized, geom="histogram",bins=100) + labs(x = 'z-score')
  ggqqplot(z_randomized) # normality of z-scores
  ks_test<-ks.test(z_randomized,'pnorm')
  p_val<-ks_test$p.value
  zscores[idx_test]=z_randomized
  print(paste("Complete",k,"of", nfolds, "folds",sep=" "))
 }

% estimate quantiles
quantiles<-c(.05, 0.25, 0.5, 0.75, 0.95)
xM<-matrix(data=NA,ncol=length(quantiles),nrow=101)
yM<-matrix(data=NA,ncol=length(quantiles),nrow=101)

xF<-matrix(data=NA,ncol=length(quantiles),nrow=101)
yF<-matrix(data=NA,ncol=length(quantiles),nrow=101)

for (i in 1:length(quantiles)){
  Qua <- getQuantile(mdl, quantile=quantiles[i],term="age",fixed.at=list(sex=0, site=300))
  c<-curve(Qua, 18, 100,  lwd=2, lty=1, add=T,col="red")
  xF[,i]<-c$x
  yF[,i]<-c$y
  
  Qua <- getQuantile(mdl, quantile=quantiles[i],term="age",fixed.at=list(sex=1, site=300))
  c<-curve(Qua, 18, 100,  lwd=2, lty=1, add=T,col="black")
  xM[,i]<-c$x
  yM[,i]<-c$y
}

fname<-"quantiles.mat"
writeMat(fname,xM=xM,yM=yM,xF=xF,yF=yF,data=data,quantiles=quantiles)
# not very good at R in plotting nice figures
# use plot_quantiles.m to generate some plots. 

 
 
