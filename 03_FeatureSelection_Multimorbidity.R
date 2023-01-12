###############################################
####  03_FeatureSelection_Multimorbidity   ####
###############################################

options(stringsAsFactors = F)

## This script performs feature selection by LASSO regression for incident "multimorbidity" 

rm(list = ls())

## set working directory 
setwd("~/rfs/rfs-epid-rfs-mpB3sSsgAn4/Studies/People/Julia/Olink_prediction_restrict_by_CVs/bin/")

## load R packages required 
library(doMC)
library(caret)
library(glmnet)
library(ROSE)

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

## generate noise variables 
rand.vars <- lapply(1:10, function(x){
  tmp <- runif(nrow(epic.p),min(sapply(epic.p[,prots], min)), max(sapply(epic.p[,prots], max)))
  return(tmp)
})

rand.vars <- as.data.frame(do.call(cbind, rand.vars))
colnames(rand.vars) <- paste0("rand.var_",1:10)

epic.sub <- cbind(epic.p, rand.vars)

## create composite "multimorbidity" outcome
epic.sub$all_disease <- rowSums(epic.sub[,d.labels.20$index_event], na.rm = T)
epic.sub$all_disease <- as.factor(ifelse(epic.sub$all_disease!=0,1,0))

## exclude prevalent cases for any of the diseases 
epic.sub$all_ex <- rowSums(epic.sub[,d.labels.20$exclusion[c(1:9,11:14)]], na.rm = T)
epic.sub$all_ex <- ifelse(epic.sub$all_ex!=0,1,0)
dat <- epic.sub[which(epic.sub$all_ex!=1),]


## ---- Run feature selection ----##

## empty objects
cf<-list()
las.morb<-list()

## draw 200 random subsets of 60% of the data
jj <- lapply(1:200,function(x) sample(nrow(dat),round(nrow(dat)*0.60),replace = F))

## kernel for parallel processing
registerDoMC(32)

## feature selection loop
for(i in 1:length(jj)){
  print(i)
  las.morb[[i]] <- train(dat[jj[[i]],c(prots, colnames(rand.vars))],
                         as.factor(dat[jj[[i]],"all_disease"]),
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
  temp<-cf[[i]][,1]
  p.select[,i]<-temp
  row.names(p.select)<-names(temp)
}

## absolute sum of weights matrix
p.select <-as.data.frame(p.select)
p.select <- p.select[-1,]
p.select$select<-rowSums(p.select[,1:ncol(p.select)])
p.select$mrc_olink.id <- row.names(p.select)
p.select <- merge(p.select, p.labels[,c("mrc_olink.id","UniProt","Assay","Panel")],by = "mrc_olink.id", all.x = T)
p.select$select.perc <- (p.select$select)/max(abs(p.select$select))
p.select$select.perc <- abs(p.select$select.perc)
    
## save feature selection ranking
write.table(p.select, file = paste0("../data_output/final_features_top20/feature_selection_all_diseases_set1_f.txt"), sep = "\t", row.names=F)

    