library(data.table) #data.table_1.13.0
library(tidyverse) #tidyverse_1.3.0
library(caret) #caret_6.0-86 
library(glmnet) #glmnet_4.0-2
library(cutpointr) #cutpointr_1.0.32
library(scales) #scales_1.1.1 
library(patchwork) #patchwork_1.0.1

#Load data needed for the function to run
data = fread("Data/Generated dataframes/LOF_all.tsv")
all_samples = fread("Data/Generated dataframes/all_donors.tsv")
features = fread("Data/COMBINED_featureCounts.tsv")
colnames(features) = str_replace_all(colnames(features), "-",".")
colnames(features) = str_replace_all(colnames(features), ">","_._")
features$Signature.1 = NULL

lasso = function(gene = "BRCA2", cohort = "Breast", Cutoff = 25){
#gene = "PBRM1"; cohort = "Kidney"; Cutoff = 25
  
  #Filter data set of mutated samples
  d = filter(data, GENE == gene)
  
  #Add mutation status for samples with no mutations in the gene
  d2 = anti_join(all_samples, d)%>%
    mutate(LOF = 0, GENE = gene, Germline = F, Double = F, MAX_CADD= 0)
  
  #Combine the two sets and filter to the cohort of interest
  d = bind_rows(d, d2)%>%
    filter( primaryTumorLocation == cohort)%>%
    dplyr::select(-cancertype)
  
  #Annotate if a sample has a pathogenic event in the gene or not
  d = d[d$MAX_CADD>Cutoff | d$MAX_CADD<5 |LOF == 7|LOF == 1 |LOF == 0, ] #Remove VUS
  d = mutate(d, y = ifelse(MAX_CADD >= Cutoff|LOF == 7, 1,0)) #Annotate pathogenic hit
  d$y[d$LOF < 2] = 0 #Annotate benign/zero hits
  
  table(d$y)
  
  #Add features
  d = left_join(d, features)
    
  #Split data into discovery (train_data) and validation (test_data)
  dfg = filter(d, y == 1)
  set.seed(1)
  nfg = sample(nrow(dfg), size = round(0.66*nrow(dfg)), replace = F)
  
  dbg = filter(d, y == 0)
  set.seed(1)
  nbg = sample(nrow(dbg), size = round(0.66*nrow(dbg)), replace = F)
  
  train_data = bind_rows(dfg[nfg,], dbg[nbg,])
  test_data = bind_rows(dfg[-nfg,], dbg[-nbg,])
  
  table(train_data$y)
  table(test_data$y)
  
  #Assign weights to the training data
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
  dplyr::select(Donor_ID, Sample_ID, LOF, GENE, y, Double, weights, Germline, LOF, primaryTumorLocation,everything())
  train_data$weights[train_data$y == 0] = w0
  
  ggplot(train_data, aes(x = weights, fill = y))+
    geom_histogram()
  
  #List of columns that needs to be removed before running the lasso, as the lasso will treat all columns as features
  cols_to_remove = c("MAX_CADD","VAF" ,"Donor_ID", "Sample_ID", "LOF","Double", "GENE", "y", "weights", "Germline", "LOF", "primaryTumorLocation")
    
  ####
  # CV.LASSO
  ####
  x = as.matrix(dplyr::select(train_data, -(cols_to_remove)))
  y = train_data$y
  
  #Set the number of cross validations — these lines make sure that each fold has both y == 1 and y == 0
  nfold <- 10
  # assign folds evenly using the modulus operator
  fold0 <- sample.int(sum(y==0)) %% nfold
  fold1 <- sample.int(sum(y==1)) %% nfold
  foldid <- numeric(length(y))
  foldid[y==0] <- fold0
  foldid[y==1] <- fold1
  foldid <- foldid + 1
  
  cv.lasso <- glmnet::cv.glmnet(x, y, 
                        alpha = 1, 
                        family = "binomial", 
                        type.measure="deviance",  
                        intercept = T,
                        foldid = foldid,
                        weights = train_data$weights, 
                        standardize = T) #We let the glmnet lasso standardize variables internally
   
  #Sanity check the binomial deviance and how it changes as lambda decreases
  plot(cv.lasso)

  #Extract the features selected by the LASSO
  f = coef(cv.lasso, s="lambda.1se")@Dimnames[[1]][coef(cv.lasso, s = "lambda.1se")@i+1] #Extract features which are kept
  
  
  # For the paper we have taken the conservative choice of using mabda = lambda.1se as generally recommended.
  # If you want a less restricted extraction of features you can however set lambda to lambda min
  #f = coef(cv.lasso, s="lambda.min")@Dimnames[[1]][coef(cv.lasso, s = "lambda.min")@i+1] #Extract features which are kept
  
  #In some cases, the LASSO returns an unreasonably large number of samples. We have capped it to a maximum of 1 feature per 10 mutated samples (rounding up to nearest 10) 
  max_f = plyr::round_any(sum(train_data$y ==1),10, ceiling) %/%10
  max_f = max_f +1 #Becasue intercept is there
  f = coef(cv.lasso, s="lambda.1se")@Dimnames[[1]][coef(cv.lasso, s = "lambda.1se")@i+1] #Extract features which are kept
  f
  if(length(f)>max_f){
    lambda_tmp = "NOPE"
    t = data.frame(nzero = cv.lasso$nzero, lambda = cv.lasso$lambda)
    while(lambda_tmp == "NOPE"){
      if(max_f %in% t$nzero){
        lambda_tmp = max(filter(t, nzero==max_f)$lambda)
      }else{
        max_f = max_f-1
      }
    }
    f = coef(cv.lasso, s=lambda_tmp)@Dimnames[[1]][coef(cv.lasso, s = lambda_tmp)@i+1] #Extract features which are kept
  }
  f = f[2:length(f)] #remove intercept
  
  #If only intercept is found, the loop ends and moves on to the next case!
  #Else, we proceed to train a GLM
  if(length(f) == 1){e_message <<- "Zero features from the lasso"; next}
  formula  = as.formula(paste("y~",paste(f, collapse = "+"), sep = ""))

  #####
  # LOOCV OF THE FINAL FORMULA ON THE TRAINING DATA
  #####
  train_data = as.data.frame(train_data)%>%
    mutate(
      y = ifelse(y == 1, "LOF", "WT")
    )
  
  test_data = as.data.frame(test_data)%>%
    mutate(
      y = ifelse(y == 1, "LOF", "WT")
    )
  
  #Now we scale and log-transform the data to be able to better interpret the output GLM
  tmp = bind_rows(mutate(train_data, s = "train"), mutate(test_data, s = "test"))
  tmp[,f] = as.numeric(scale(log(tmp[,f]+0.1))) #Scale it allll together :D
  for(col in f){ 
    tmp_col = dplyr::select(tmp, col)
    tmp_col = rescale(tmp_col[,col], to = c(0,1)) #Rescale to get values pushed to the positive interval 0—1
    tmp[,col] = tmp_col
  }
  
  #Change the "y" to a string of either LOF (loss-of-function) or WT (wildtype)
  tmp$y= factor(tmp$y, levels = c("WT", "LOF")) #Set the WT as the 'positive' class to get estimates for the transition into WT
  
  train_data = filter(tmp, s == "train")%>%
    dplyr::select(-s)
  test_data = filter(tmp, s != "train")%>%
    dplyr::select(-s)
  
  #Train a GLM on the training data
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
  
  #Some models have very high std.errors, as the GLM fails to properly converge. In these cases, we jump to the next case
  if(max(t$`Std. Error`) >= 20){next}
  
  #Else, we plot the coefficients
  Features = data.frame(Gene = gene, Cohort = cohort, 
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
  
  b1 = table(tmp$y)
  
  p1 = ggplot(Features, aes(x= Estimate, y = Feature))+
    geom_errorbar(aes(xmin = Estimate - stderror, xmax = Estimate + stderror))+
    theme_bw(base_size = 9)+
    theme(
     # legend.position = c(0.75, 0.1),
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
  train_set = dplyr::select(train_set, -weights, -rowIndex, -parameter)
  
  test_set = predict.train(c, test_data, type = "prob")
  test_set$"obs" = factor(test_data$y, levels = c("WT", "LOF")) #Set the WT as the 'positive' class to get estimates for the transition into WT
  test_set$"pred" = as.factor(predict.train(c, test_data, type = "raw"))
  test_set$mutation_level = test_data$LOF
  test_set = test_set%>%
    mutate(
      Source = "Testing"
    )
  
  #########
  ####Evaluate the Predictive power
  #########
  cp_train = cutpointr(train_set, x = LOF, class = obs, pos_class = "LOF",
                       method = maximize_metric, metric = sum_sens_spec,  direction = ">=")
  train = as.data.frame(cp_train$roc_curve)
  train$"label" = paste("Discovery, AUC: ", round(cp_train$AUC,3), sep ="")
  
  genenames = paste(gene, collapse = "_")
  cohortnames = paste(cohort, collapse = "_")
  
  #Here you may save the model, with the optimal cutpoint in the filename
  #saveRDS(object = model_out, file = paste("Results/Models/",genenames,cohortnames,cp_train$optimal_cutpoint,".Rdata" ,sep = "|"))
  
  cp_test = cutpointr(test_set, x = LOF, class = obs, pos_class = "LOF",
                      method = oc_manual, cutpoint = cp_train$optimal_cutpoint, direction = ">=")
  test = as.data.frame(cp_test$roc_curve)
  test$"label" = paste("Validation, AUC: ", round(cp_test$AUC,3), sep ="")
  
  #Make ROC plot  
  tt = bind_rows(train, test)
  p2 = ggplot()+
    geom_line(data = tt, aes(x = 1-tnr, y = tpr, color = label))+
    theme_bw(base_size = 9)+
    theme(
     # legend.position = c(0.8, 0.2),
      legend.position = "bottom",
      legend.direction = "vertical",
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
  ggsave(plot = p3, paste("Results/Models/",gene, cohort,".pdf", sep = ""), device = "pdf", width = 7, height = 3)
  

  #####
  # Combine and return results
  #####

  full_model=
    sprintf('y ~ %.2f + %s', coef(c$finalModel)[1], 
            paste(sprintf(' (%.2f)*%s ',
                          coef(c$finalModel)[-1], names(coef(c$finalModel)[-1]) ), collapse ='+'))
  
  Model = data.frame(
    Gene = paste(gene, collapse = "|"),
    Cohort = paste(cohort, collapse = "|"),
    Train_fg = nrow(train_data[train_data$y == "LOF",]),
    Train_bg = nrow(train_data[train_data$y == "WT",]),
    Test_fg = nrow(test_data[test_data$y == "LOF",]),
    Test_bg = nrow(test_data[test_data$y == "WT",]),
    Formula = full_model,
    ROC = cp_train$AUC,
    ROC_test = cp_test$AUC,
    Sensitivity = cp_train$sensitivity,
    Specificity = cp_train$specificity
  )

  
  print("Done Model")
  
  output = list(Model, Features)
  return(output)
}

