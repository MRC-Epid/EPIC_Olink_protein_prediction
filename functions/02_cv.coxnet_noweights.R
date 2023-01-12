## function to run coxnet and get cross-validation C-index
coxnet.cv <- function(Train.data, predictors, Train.surv.data,times ) {
  ##  Train.data : data set used for training
  ##  predictors : vector of predictors
  ##  Train.surv.data : survival obset from the Training data
  ##  cox.weights : optional weigths to be used
  ##  times : number of boots for getting the cindex stats
  ##  Test.data : data set used for Testing
  ##  Test.surv.data : survival objects in the testing data
  
  Train.surv.data <- Train.surv.data[complete.cases(Train.data[,predictors]),]
  Train.data <- Train.data[complete.cases(Train.data[,predictors]),]
  
  
  require(caret)
  require(glmnet)
  require(ROSE)
  ## set seed and optimize
  boot.cidenx <- NULL
  lo.cidenx <- NULL
  up.cidenx <- NULL
  las.morb <- list()
  for(i in 1:times){
    las.morb[[i]] <- cv.glmnet(as.matrix(Train.data[,predictors]),
                          as.matrix(Train.surv.data),
                          family="cox",alpha=0,type.measure = "C",
                          lambda=10^-seq(10,.25,-.25),
                          nfolds = 5,sampling="rose"
    )
    index <- las.morb[[i]]$index[1]
    boot.cidenx <- c(boot.cidenx, las.morb[[i]]$cvm[index])
    lo.cidenx <- c(lo.cidenx, las.morb[[i]]$cvlo[index])
    up.cidenx <- c(up.cidenx, las.morb[[i]]$cvup[index])
  }
  glmnet.opt <- las.morb[[which.max(boot.cidenx)]]$glmnet.fit
  s.opt <- las.morb[[which.max(boot.cidenx)]]$lambda.min
  opt.coef <- coef(glmnet.opt,s = s.opt)
  coxnet.res <- list(glmnet.opt=las.morb[[which.max(boot.cidenx)]]$glmnet.fit,
                     lambda.opt=s.opt,coefficients=opt.coef,
                     Cindex.vec=boot.cidenx, Cindex.mean=mean(boot.cidenx),
                     ci.low=mean(lo.cidenx), ci.upper=mean(up.cidenx))
  return(coxnet.res)
}
  