library(data.table)
library(tidyverse)
library(caret)
library(glmnet)
library(patchwork)
library(scales)
library(extrafont)
#library(gridExtra)
library(cutpointr)
setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/R-projects/DDR_project_Jan_2021/")

data = fread("Data/LOF_all.tsv")
all_samples = fread("Data/all_donors.tsv")
features = fread("Data/COMBINED_featureCounts.tsv")
colnames(features) = str_replace_all(colnames(features), "-",".")
colnames(features) = str_replace_all(colnames(features), ">","_._")
data = filter(data, Sample_ID %in% features$Sample_ID)
all_samples = filter(all_samples, Sample_ID %in% features$Sample_ID)
table(!duplicated(data$Sample_ID))
features$Signature.1 = NULL

  gene = "BRCA1"; cohort = c("Ovary", "Breast"); Cutoff = 25
  
  #Filter data set
  d = filter(data, GENE == "BRCA1")
  d2 = anti_join(all_samples, d)%>%
    mutate(LOF = 0, GENE = gene, Germline = F, Double = F, MAX_CADD= 0)
  d = bind_rows(d, d2)%>%
    filter( primaryTumorLocation %in% cohort)%>%
    dplyr::select(-cancertype)
  d = d[d$MAX_CADD>Cutoff | d$MAX_CADD<5 |LOF == 7|LOF == 1 |LOF == 0, ]
  d = mutate(d, y = ifelse(MAX_CADD >= Cutoff|LOF == 7, 1,0))
  d$y[d$LOF < 2] = 0
  
  #Remove BRCA2
  brca2 = filter(data, GENE == "BRCA2")
  brca2 = filter(brca2, LOF >= 6)
  
  d = filter(d, !Sample_ID %in% brca2$Sample_ID)
  
  table(d$y)
  
  #Add features
  d = left_join(d, features)
  
  train_data = d
  
  #  table(train_data$y)
  # table(test_data$y)
  
  #Assign weights to the training
  w = table(train_data$y)
  w1 = w[1]/nrow(train_data)
  w0 = w[2]/nrow(train_data)
  train_data = mutate(train_data,
                      y = as.factor(y),
                      weights = ifelse(y == 1, w1, w0),
                      weights = weights + w1* (MAX_CADD-30)/50  , #E.g. CADD 32 = 32-30 = 2, 2/20, 0.1, 0.1*WD            
                      weights = ifelse(!is.na(VAF), weights + w1* (VAF-35)/100, weights),   #E.g. VAF == 40 = 40-35 = 5, 5/100 = 0.05*WD                      weights = ifelse(Double == T, (weights * 1.20), weights),
                      weights = ifelse(LOF == 7, weights + w1 * 0.20, weights), #If ClinVar add 20%
                      weights = ifelse(Double == T, weights + w1 * 0.20, weights) #If double hit add 20%
  )%>% 
    dplyr::select(Donor_ID, Sample_ID, LOF, GENE, y, Double, weights, Germline, LOF, Study, primaryTumorLocation,everything())
  
  train_data$weights[train_data$y == 0] = w0
  
  hist(train_data$weights, breaks = 1e2)
  
  cols_to_remove = c("MAX_CADD","Donor_ID", "VAF", "Sample_ID", "LOF","Double", "GENE", "y", "weights", "Germline", "LOF", "Study", "primaryTumorLocation")
  
  ####
  # CV.LASSO
  ####
  x = as.matrix(dplyr::select(train_data, -(cols_to_remove)))
  y = train_data$y
  
  nfold <- 30
  # assign folds evenly using the modulus operator
  fold0 <- sample.int(sum(y==0)) %% nfold
  fold1 <- sample.int(sum(y==1)) %% nfold
  foldid <- numeric(length(y))
  foldid[y==0] <- fold0
  foldid[y==1] <- fold1
  foldid <- foldid + 1
  
  #  tryCatch({
  cv.lasso <- cv.glmnet(x, y, 
                        alpha = 1, 
                        family = "binomial", 
                        type.measure="deviance",  
                        intercept = T,
                        foldid = foldid,
                        #nfolds = nrow(x),
                        weights = train_data$weights,
                        standardize = T) #Use CV Lasso to estimate best model and return it
  # }, warning=function(cond){e_message <<- "Lasso doesn't converge";next})
  
  plot(cv.lasso)
  
  #Max_features
  max_f = plyr::round_any(sum(train_data$y ==1),10, ceiling) %/%10
  max_f = max_f +1
  f = coef(cv.lasso, s="lambda.min")@Dimnames[[1]][coef(cv.lasso, s = "lambda.min")@i+1] #Extract features which are kept
  # f
  # if(length(f)>max_f){
  #   lambda_tmp = "NOPE"
  #   t = data.frame(nzero = cv.lasso$nzero, lambda = cv.lasso$lambda)
  #   while(lambda_tmp == "NOPE"){
  #     if(max_f %in% t$nzero){
  #       lambda_tmp = max(filter(t, nzero==max_f)$lambda)
  #     }else{
  #       max_f = max_f-1
  #     }
  #   }
  #   f = coef(cv.lasso, s=lambda_tmp)@Dimnames[[1]][coef(cv.lasso, s = lambda_tmp)@i+1] #Extract features which are kept
  # }
  f = f[2:length(f)]
  f
  cat("N_features:", length(f))
  
  if(length(f) < 1){e_message <<- "Zero features from the lasso"; next}
  
  print("Time to validate!")
  
  formula  = as.formula(paste("y~",paste(f, collapse = "+"), sep = ""))
  
  #####
  # LOOCV OF THE FINAL FORMULA ON THE TRAINING DATA
  #####
  train_data = as.data.frame(train_data)%>%
    mutate(
      y = ifelse(y == 1, "LOF", "WT")
    )
  
  train_data[,f] = as.numeric(scale(log(train_data[,f]+0.1))) #Scale it allll together :D
  for(col in f){ 
    tmp_col = dplyr::select(train_data, col)
    tmp_col = rescale(tmp_col[,col], to = c(0,1)) #Rescale to get values pushed to the positive interval 0—1
    train_data[,col] = tmp_col
  }
  
  train_data$y= factor(train_data$y, levels = c("WT", "LOF")) #Set the WT as the 'positive' class to get estimates for the transition into WT
  
  c = caret::train(formula, method = "glm", 
                   data = train_data, 
                   trControl = trainControl(method = "LOOCV", classProbs = TRUE),
                   #    metric = "RMSE",
                   weights = weights,
                   family = "quasibinomial"
  )
  
  summary(c)
  model_out = c
  
  t = as.data.frame(summary(c)$coefficients)
  Features = data.frame(Gene = gene, Cohort = paste(cohort, collapse = ", "), 
                        Feature = rownames(t), 
                        Estimate = t$Estimate,
                        stderror = t$`Std. Error`,
                        t_val = t$`t value`, 
                        p_val = t$`Pr(>|t|)`)
  Features$Feature = str_replace_all(Features$Feature, "\\__._1", ">1")
  Features$Feature = str_replace_all(Features$Feature, "b\\.", "b-")
  Features$Feature = str_replace_all(Features$Feature, "1.1", "1-1")
  Features$Feature = str_replace_all(Features$Feature, "\\.", " ")
  Features$Feature = str_replace_all(Features$Feature, "del", "del.")
  Features$Feature = str_replace_all(Features$Feature, "rep", "rep.")
  Features$Feature = str_replace_all(Features$Feature, "mh", "Microhomology")
  Features$Feature = str_replace_all(Features$Feature, "\\_", " ")
  
  b1 = table(train_data$y)
  
  p1 = ggplot(Features, aes(x= Estimate, y = Feature))+
    geom_errorbar(aes(xmin = Estimate - stderror, xmax = Estimate + stderror))+
    theme_bw(base_size = 12, base_family = "Arial")+
    theme(
      legend.position = c(0.75, 0.1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      axis.line = element_line(color = NA)
    )+
    geom_vline(xintercept = 0, lty = 2, col = "grey")+
    ylab("")+
    ggtitle(paste(gene, ", ", cohort, " cancer",sep = ""), 
            subtitle = paste( "WT: ", b1[1], ", LOF: ", b1[2], " (Train: ",sum(train_data$y == "LOF"), ")", sep =""))
  p1
  
  #####
  # Predict on testing data
  #####
  train_set = c$pred
  train_set$mutation_level = train_data$LOF 
  train_set$"Source" = "Training"
  # train_set$"Study" = study
  train_set = dplyr::select(train_set, -weights, -rowIndex, -parameter)
  train_set$CT = train_data$primaryTumorLocation
  
  #########
  ####Predictive power
  #########
  cp_train = cutpointr(train_set, x = LOF, class = obs, pos_class = "LOF",
                       method = maximize_metric, metric = sum_sens_spec,  direction = ">=", subgroup = CT)
  
  for(j in 1:length(cohort)){
    tmpRoc = as.data.frame(cp_train$roc_curve[[j]])
    tmpRoc$"label" = paste(cohort[[j]], ", AUC: ", round(cp_train$AUC[[j]],3), sep ="")
    if(j == 1){
      train = tmpRoc
    }else{
      train = bind_rows(train, tmpRoc)
    }
  }
  
  genenames = paste(gene, collapse = "_")
  cohortnames = paste(cohort, collapse = "_")
  cutoffs =  paste(cp_train$optimal_cutpoint, collapse = "_")
  saveRDS(object = model_out, file = paste("Results/Valid_Models/",genenames,cohortnames,cutoffs,".Rdata" ,sep = "|"))
  
  #Make ROC plot  
  p2 = ggplot()+
    geom_line(data = train, aes(x = 1-tnr, y = tpr, color = label))+
    theme_bw(base_size = 12, base_family = "Arial")+
    theme(
      legend.position = c(0.8, 0.2),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      axis.line = element_line(color = NA)
    )+
    xlab("1-Specificity")+
    ylab("Sensitivity")+
    scale_color_manual(values = c("Darkblue", "brown"))+
    guides(color = guide_legend(title = ""))
  p2
  
  p3=p1+p2 
  p3
  ggsave(plot = p3, paste("Results/ROCS_Valid/",gene, cohort,".pdf", sep = ""), device = "pdf", width = 13, height = 4)
  
  
  #####
  # Combine and return results
  #####
  print("Make Model")
  
  full_model=
    sprintf('y ~ %.2f + %s', coef(c$finalModel)[1], 
            paste(sprintf(' (%.2f)*%s ',
                          coef(c$finalModel)[-1], names(coef(c$finalModel)[-1]) ), collapse ='+'))
  
  # conf = as.data.frame(conf$table)
  # conf$Gene = gene
  # conf$cohort = cohort
  
  Model = data.frame(
    Gene = paste(gene, collapse = "|"),
    Cohort = paste(cohort, collapse = "|"),
    # Training_study = study,
    Train_fg = nrow(train_data[train_data$y == "LOF",]),
    Train_bg = nrow(train_data[train_data$y == "WT",]),
    Formula = full_model,
    ROC = paste(cp_train$AUC, collapse = "|"),
    Sensitivity = paste(cp_train$sensitivity, collapse = "|"),
    Specificity = paste(cp_train$specificity, collapse = "|"),
    p_val_study =   wilcox.test(x = train_set$LOF[train_set$obs == "LOF"], 
                                y = train_set$LOF[train_set$obs == "WT"], 
                                alternative = "greater", exact = F)$p.value,
    p_val_test = NA
  )
  
  
  print("Done Model")
  
  output = list(Model, Features)
  # print("Returning output")
  return(output)
}

#lasso(gene = gene, germline = F, cohort = cohort, study = study, remove_samples = "NO", collectSV = F)
#lasso()[[3]]

#o = lasso(gene = "BRCA2", cohort = "Breast")
#o
#rm(d_lof)
#out = lasso(reps = 100, gene = "BRCA2", germline = TRUE)  


