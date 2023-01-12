##############################################
####  04_CrossValidation_Multimorbidity   ####
##############################################

options(stringsAsFactors = F)

## This script runs validation of predictive multimorbidity protein signature for individual diseases by cross-validation   

rm(list = ls())

## set working directory 
setwd("~/rfs/rfs-epid-rfs-mpB3sSsgAn4/Studies/People/Julia/Olink_prediction_restrict_by_CVs/bin/")

## load R packages required
library(doMC)
library(caret)
library(glmnet)
library(ROSE)
library(dplyr)
library(tidyr)
library(survival)

## read in cleaned batch 2 data
epic.p <- read.delim("../../TATAA_validation_fixONCII/data_input/EPIC-tataa3K-phenos_subcohort.txt", sep = "\t")

## read in Olink protein metadata
p.labels <- read.delim("../../Olink_explore_TATAA_replication/data_input/TATAA_3K_EPICN_labels.txt", sep = "\t")
p.labels$mrc_olink.id <- gsub("-",".", p.labels$mrc_olink.id)

## restrict to proteins with comparable coefficients of variantion across the two batches
cvs <- read.delim("../data_input/cv_set1_set_2.txt", sep = "\t")
cvs.stable <- subset(cvs, color==0)
prots <- p.labels$mrc_olink.id[which(p.labels$OlinkID%in%cvs.stable$OlinkID)]

## read in disease labels prevalent disease labels and exclusion criteria labels
d.labels.20 <- read.delim("../../TATAA_validation_fixONCII/data_input/disease_labels_set2_subcohort.txt", sep = "\t")

## PRS name list per disease 
grs.list <- list("inc_diabetes"="grs.t2d_s","inc_renadis"="grs.ckd_s","inc_ihd"="grs.cad_s","inc_af"="grs.af_s","inc_cardfail"="grs.cardfail_s",
                 "inc_cvainfct"="grs.cvainfct_s","inc_pad_all"="grs.pad_s","inc_venous_t"="grs.venous_t_s","inc_resp_asth"="grs.asthma_s",
                 "inc_resp_copd"="grs.copd_s","breast_ca"="grs.breast_ca_s","nmskin_ca"="grs.skin_ca_s",
                 "dead"="grs.telomere_s","prostate_ca"="grs.prostate_ca_s","cataract_hosp"="grs.cataracts_s","glaucoma_hosp"="grs.glaucoma_s",
                 "frac_all_hosp"="grs.fractures_s")

epic.sub <- epic.p

## read in optimization function
source("functions/02_cv.coxnet_noweights.R")

## Run cross-validation across diseases
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
  
  ## recode to compute premature mortality (mortality before age of 75) 
  if(d.labels.20$index_event[x]=="dead"){
    dat$age_stop_FU <- dat$age+dat$fol
    
    dat <- dat[which((dat[d.labels.20$index_event[x]]==1&dat$age_stop_FU<75)|dat[d.labels.20$index_event[x]]==0),]
    
  }
  
  ## generate survival object
  surv.dat <- Surv(dat$fol,dat[,d.labels.20$index_event[x]])
  
  ## kernel for parallel processing
  registerDoMC(32)
  
  ## print disease currently running 
  print(d.labels.20$index_event[x])
  
  ## Read in top 20 multimorbidity protein signature
  fs.20 <- read.delim(paste0("../data_output/final_features_top20/feature_selection_all_diseases_set1_f.txt"), sep = "\t")
  
  ## restrict to top 10
  predictors <- fs.20 %>%
    top_n(10,select.perc)
  predictors <- as.character(predictors$mrc_olink.id)
  
  ## define patient information predictor variables
  if(!d.labels.20$index_event[x]%in%c("breast_ca","prostate_ca")){
    clin.vars <- c("age","sex","bmi","ever_smoker")
  } else{
    clin.vars <- c("age","bmi","ever_smoker")
  }
  
  ## run cross-validation over bootstraping
  protein.coxnet <- coxnet.cv(Train.data = dat,
                                predictors =predictors ,
                                Train.surv.data = surv.dat,
                                times =100 )
  clinical.coxnet <- coxnet.cv(Train.data = dat,
                                 predictors = clin.vars ,
                                 Train.surv.data = surv.dat,
                                 times =100 )
  clinical.protein.coxnet <- coxnet.cv(Train.data = dat,
                                       predictors =c(predictors, clin.vars) ,
                                       Train.surv.data = surv.dat,
                                       times =100 )

  ## save results in a list 
  res.list <- list(protein.model=protein.coxnet,
                   clin.model=clinical.coxnet,
                   clin.protein.model = clinical.protein.coxnet,
                   n=nrow(dat),
                   n.cases=sum(dat[,d.labels.20$index_event[x]]==1, na.rm = T)
  )
  return(res.list)
})
## Disease names to list 
names(lasso) <- d.labels.20$index_event

## extract average c-index results
cindex <- t(sapply(d.labels.20$index_event, function(x){
  return(data.frame( protein.model= lasso[[x]]$protein.model$Cindex.mean,
                     clinical.model= lasso[[x]]$clin.model$Cindex.mean,
                     clinical.protein.model= lasso[[x]]$clin.protein.model$Cindex.mean ))
}))
cindex <- as.data.frame(cindex)
cindex$disease <- row.names(cindex)

## extract average lower error curve
low.ci <- t(sapply(d.labels.20$index_event, function(x){
  return(data.frame( protein.model= lasso[[x]]$protein.model$ci.low,
                     clinical.model= lasso[[x]]$clin.model$ci.low,
                     clinical.protein.model= lasso[[x]]$clin.protein.model$ci.low ))
}))
low.ci <- as.data.frame(low.ci)
low.ci$disease <- row.names(low.ci)

## extract average upper error curve
up.ci <- t(sapply(d.labels.20$index_event, function(x){
  return(data.frame( protein.model= lasso[[x]]$protein.model$ci.upper,
                     clinical.model= lasso[[x]]$clin.model$ci.upper,
                     clinical.protein.model= lasso[[x]]$clin.protein.model$ci.upper ))
}))
up.ci <- as.data.frame(up.ci)
up.ci$disease <- row.names(up.ci)

## long format to merge
cindex <- gather(cindex, model, mean.Cindex, 1:3)
low.ci <- gather(low.ci, model , low.ci, 1:3)
up.ci <- gather(up.ci, model , up.ci, 1:3)

## merge all
ci.res <- merge(cindex, low.ci, by = c("disease","model"))
ci.res <- merge(ci.res, up.ci, by = c("disease","model"))

## to numeric
ci.res$mean.Cindex <- as.numeric(ci.res$mean.Cindex)
ci.res$low.ci <- as.numeric(ci.res$low.ci)
ci.res$up.ci <- as.numeric(ci.res$up.ci)

## save results table
write.table(ci.res, file = "../data_output/cindex_FS_all_diseases_top10_sc1_cvopti_sc2.txt", sep = "\t", row.names = F)
