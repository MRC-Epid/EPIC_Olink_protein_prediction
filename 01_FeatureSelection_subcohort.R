###########################################
####  01_Feature_Selection_subcohort   ####
###########################################

options(stringsAsFactors = F)

## This script performs feature selection by LASSO regression for incident disease status    

rm(list = ls())

## set working directory 
setwd("~/rfs/rfs-epid-rfs-mpB3sSsgAn4/Studies/People/Julia/Olink_prediction_restrict_by_CVs/bin/")

## load R packages required 
library(doMC)
library(caret)
library(glmnet)
library(ROSE)
library(survival)
library(dplyr)
library(tidyr)

## read in cleaned batch 1 data
epic.p <- read.delim("../../Olink_diseases_clean_de_novo/data_input/EPIC-olink3K-phenos_subcohort.txt", sep = "\t")

## read in Olink protein metadata
p.labels <- read.delim("../../Olink_explore_TATAA_replication/data_input/TATAA_3K_EPICN_labels.txt", sep = "\t")
p.labels$mrc_olink.id <- gsub("-",".", p.labels$mrc_olink.id)

## restrict to proteins with comparable coefficients of variantion across the two batches
cvs <- read.delim("../data_input/cv_set1_set_2.txt", sep = "\t")
cvs.stable <- subset(cvs, color==0)
prots <- p.labels$mrc_olink.id[which(p.labels$OlinkID%in%cvs.stable$OlinkID)]

## read in disease labels prevalent disease labels and exclusion criteria labels
d.labels.20 <- read.delim("../../Olink_diseases_clean_de_novo/data_input/disease_labels_set1_subcohort.txt", sep = "\t")

## PRS name list per disease
grs.list <- list("inc_diabetes"="grs.t2d_s","inc_renadis"="grs.ckd_s","inc_ihd"="grs.cad_s","inc_af"="grs.af_s","inc_cardfail"="grs.cardfail_s",
                 "inc_cvainfct"="grs.cvainfct_s","inc_pad_all"="grs.pad_s","inc_venous_t"="grs.venous_t_s","inc_resp_asth"="grs.asthma_s",
                 "inc_resp_copd"="grs.copd_s","breast_ca"="grs.breast_ca_s","nmskin_ca"="grs.skin_ca_s",
                 "dead"="grs.telomere_s","prostate_ca"="grs.prostate_ca_s","cataract_hosp"="grs.cataracts_s","glaucoma_hosp"="grs.glaucoma_s",
                 "frac_all_hosp"="grs.fractures_s")

## Generate random noise variables 
rand.vars <- lapply(1:10, function(x){
  tmp <- runif(nrow(epic.p),min(sapply(epic.p[,prots], min)), max(sapply(epic.p[,prots], max)))
  return(tmp)
})

rand.vars <- as.data.frame(do.call(cbind, rand.vars))
colnames(rand.vars) <- paste0("rand.var_",1:10)

epic.sub <- cbind(epic.p, rand.vars)

## read in optimization function
source("functions/01_coxnet_ridge_optim_noweights.R")

lasso <- lapply(1:nrow(d.labels.20), function(x){
  ## exclusion of prevalent disease cases
  ## for sex specific outcomes, exclusion of other sex
  if (d.labels.20$index_event[x]%in%c("inc_diabetes","inc_renadis","inc_ihd",
                                      "inc_af","inc_cardfail","inc_cvainfct",
                                      "inc_pad_all","inc_venous_t","inc_resp_asth",
                                      "inc_resp_copd","nmskin_ca","cataract_hosp",
                                      "glaucoma_hosp","frac_all_hosp" )) {
    dat <- epic.sub[which(epic.sub[d.labels.20$exclusion[x]]!=1),]
  }
  if (d.labels.20$index_event[x]=="dead") {
    dat <- epic.sub
  }
  if (d.labels.20$index_event[x]=="prostate_ca") {
    dat <- epic.sub[which(epic.sub[d.labels.20$exclusion[x]]!=1&epic.sub$sex==1),]
  }
  if (d.labels.20$index_event[x]=="breast_ca") {
    dat <- epic.sub[which(epic.sub[d.labels.20$exclusion[x]]!=1&epic.sub$sex==2),]
  }
  
  ## compute follow-up time
  dat$fol = dat[,d.labels.20$index_date[x]] - dat[,"ndate"]
  
  ## subset to incident cases after 6 months
  dat <- subset(dat, fol>=0.5)
  dat <- dat[complete.cases(dat[,d.labels.20$index_event[x]]),]
  
  ## generate survival object
  surv.dat <- Surv(dat$fol,dat[,d.labels.20$index_event[x]])
  
  ## ---- Run feature selection ----##
  
  ## empty objects
  cf<-list()
  las.morb<-list()

  ## draw 200 random subsets of 60% of the data
  jj <- lapply(1:200,function(x) sample(nrow(dat),round(nrow(dat)*0.60),replace = F))

  ## kernel for parallel processing
  # registerDoMC(32)
  
  ## print disease running 
  print(d.labels.20$index_event[x])
  ## Feature selection loop 
    for(i in 1:length(jj)
        ){
        print(i)
        las.morb[[i]] <- train(dat[jj[[i]],c(prots, colnames(rand.vars))],
                               as.factor(dat[jj[[i]],d.labels.20$index_event[x]]),
                               metric = "Kappa",
                               family="binomial",
                               method = "glmnet",
                               tuneGrid = as.data.frame(expand.grid(alpha=1, lambda=10^-seq(10,.25,-.25))),
                               trControl = trainControl(method="repeatedcv", number=3, repeats = 5,
                                                        sampling="rose", allowParallel=T)
        )
        cf[[i]]       <- as.matrix(coef(las.morb[[i]]$finalModel, s = las.morb[[i]]$finalModel$lambdaOpt))
        las.morb[[i]] <- las.morb[[i]]$results
        gc()
    }
   
    ## Feature selection weights matrix
    p.select <- matrix(, nrow = nrow(cf[[1]]),ncol = length(cf))
    
    for (i in 1:length(cf)){
      temp <- cf[[i]][,1]
      p.select[,i] <- temp
      row.names(p.select) <- names(temp)
    }
    
    ## absolute sum of weights matrix
    p.select <-as.data.frame(p.select)
    p.select <- p.select[-1,]
    p.select$select<-rowSums(p.select[,1:ncol(p.select)])
    p.select$mrc_olink.id <- row.names(p.select)
    p.select <- merge(p.select, p.labels[,c("mrc_olink.id","UniProt","Assay","Panel")],by = "mrc_olink.id", all.x = T)
    p.select$select.perc <- (p.select$select)/max(abs(p.select$select))
    p.select$select.perc <- abs(p.select$select.perc)
    
    ## check highest selection score for random variables
    t.sel <- p.select$select.perc[grep("rand", p.select$mrc_olink.id)[1]]
    
    ## select top 20 predictors
    predictors <- p.select %>%
      top_n(20, select.perc)
    predictors <- predictors$mrc_olink.id
    
    ## define clinical variables
    if(!d.labels.20$index_event[x]%in%c("breast_ca","prostate_ca")){
      clin.vars <- c("age","sex","bmi","ever_smoker")
    } else{ 
      clin.vars <- c("age","bmi","ever_smoker")
      }
    
    ## optimization using top 20 proteins 
    sel.weigths <- coxnet.optim.r(Train.data = dat,predictors =c(predictors,clin.vars) ,
                                Train.surv.data = surv.dat,
                                times =1 , Test.data =dat , Test.surv.data = surv.dat)
    
    ## weights from reduced model includeing top 20 proteins and clinical variables
    predictors <- sel.weigths$coefficients[grep("_[0-9]{5}", row.names(sel.weigths$coefficients)),]
    predictors <- predictors %>%
      as_tibble(rownames = "mrc_olink.id") %>%
      mutate(abs.coef = abs(value)) %>%
      filter(abs.coef!=0)
    predictors <- merge(predictors, p.select[,c("mrc_olink.id","select","select.perc")], by = "mrc_olink.id")
    ## rank by final selection score (product of feature selection and optimization weights)
    predictors <- predictors%>%
      mutate(sel.opti.score = abs.coef*abs(select))
    
    write.table(predictors, file = paste0("../data_output/final_features_top20/",d.labels.20$index_event[x],"_set1.txt"), sep = "\t", row.names=F)
})
    