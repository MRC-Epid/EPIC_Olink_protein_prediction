## function to run optimization by LASSO CV ## normal cox without PW-weights
coxnet.optim.r <- function(Train.data, predictors, Train.surv.data,times,Test.data,Test.surv.data ) {
  ##  Train.data : data set used for training
  ##  predictors : vector of predictors
  ##  Train.surv.data : survival obset from the Training data
  ##  cox.weights : optional weigths to be used
  ##  times : number of boots for getting the cindex stats
  ##  Test.data : data set used for Testing
  ##  Test.surv.data : survival objects in the testing data
  
  Train.surv.data <- Train.surv.data[complete.cases(Train.data[,predictors]),]
  Test.surv.data <- Test.surv.data[complete.cases(Test.data[,predictors]),]
  Train.data <- Train.data[complete.cases(Train.data[,predictors]),]
  Test.data <- Test.data[complete.cases(Test.data[,predictors]),]
  
  
  require(caret)
  require(glmnet)
  require(ROSE)
  ## set seed and optimize
  set.seed(354)
  las.morb <- cv.glmnet(as.matrix(Train.data[,predictors]),
                             as.matrix(Train.surv.data),
                             family="cox",alpha=0,type.measure = "C",
                             lambda=10^-seq(10,.25,-.25),
                             nfolds = 5 #,sampling="rose"
                        )
  ## save the model and parameters from best CV round
  glmnet.opt <- las.morb$glmnet.fit
  s.opt <- las.morb$lambda.min
  opt.coef <- coef(glmnet.opt,s = s.opt)
  ## Generate Cindex statostics by bootstrapping
  boot.cidenx <- NULL
  jj <- list()
  tmp.pred <- list()
  tmp.index <- list()
  for (i in 1:times){
    jj[[i]] <- sample(nrow(Test.data),nrow(Test.data),replace = T)
    tmp.pred[[i]] <- as.numeric(predict(glmnet.opt,type = "response",newx = as.matrix(Test.data[jj[[i]],predictors]),s=s.opt))
    tmp.index[[i]] <- glmnet::Cindex(tmp.pred[[i]],as.matrix(Test.surv.data[jj[[i]]]))
    boot.cidenx <- c(boot.cidenx,tmp.index[[i]])
  }
  ## aggregate all results in a list
  coxnet.res <- list(glmnet.opt=glmnet.opt,lambda.opt=s.opt,coefficients=opt.coef,Cindex.vec=boot.cidenx, Cindex.mean=mean(boot.cidenx),
                     ci.low=quantile(boot.cidenx,0.025), ci.upper=quantile(boot.cidenx,0.975))
  return(coxnet.res)
}
