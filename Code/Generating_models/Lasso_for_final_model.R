
lassoLoop = function(train_data = train_data, test_data = test_data, nIterations = 10, nfold = 5, Null = F, 
                     test = F,cols_to_remove = cols_to_remove, gene = gene, cohort = cohort, study = trainStudy){
  
  #train_data = train_data; test_data = test_data; nIterations = 10; nfold = 5; Null = F; test = F;cols_to_remove = cols_to_remove; gene = gene; cohort = cohort; study = trainStudy
  for(i in 1:nIterations){
  
  if(!cols_to_remove == "None"){
    x = dplyr::select(train_data, -c(cols_to_remove))
  }
  x = as.matrix(x)
  y = train_data$y
  weights = train_data$weights
  
  if(Null == T){ #Reorder which samples are "mutated" and make sure they are equally weighted
    tmp_y = data.frame(y = y, weights = weights) 
    tmp_y = tmp_y[sample(1:nrow(tmp_y), size = nrow(tmp_y), replace =F),]
    weights =tmp_y$weights
    y = tmp_y$y
    #print(table(y))
  }
  
  set.seed(i)
  nfold <- nfold
  # assign folds evenly using the modulus operator
  fold0 <- sample.int(sum(y==0)) %% nfold
  fold1 <- sample.int(sum(y==1)) %% nfold
  foldid <- numeric(length(y))
  foldid[y==0] <- fold0
  foldid[y==1] <- fold1
  foldid <- foldid + 1
  
  #  tryCatch({
  cv.lasso <<- glmnet::cv.glmnet(
    x = x, 
    y = factor(y), 
    alpha = 1, 
    family = "binomial", 
    type.measure="deviance",   
    intercept = T,
    foldid = foldid,
    #nfolds = nrow(x), #Kfold cross validation, performs worse!
    weights = train_data$weights, 
    #lower.limits= 0,
    standardize = F) #Use CV Lasso to estimate best model and return it
  # }, warning=function(cond){e_message <<- "Lasso doesn't converge";next})
  #plot(cv.lasso)
  #We use lambda = 1SE
 # summary(cv.lasso$cvm) #See AUC scores of the runs
  
  #Cap the set of features
  max_f = plyr::round_any(sum(train_data$y ==1),10, ceiling) %/%10
  max_f = max_f + 1 #+1 for the intercept
  if(max_f > 4){max_f = 4}
  #max_f = 6 #Max 5 features next to the intercept
  f = coef(cv.lasso, s="lambda.min")@Dimnames[[1]][coef(cv.lasso, s = "lambda.min")@i+1] #Extract features which are kept
  f
  if(length(f)>max_f){
    lambda_tmp = "NOPE"
    t = data.frame(nzero = cv.lasso$nzero, lambda = cv.lasso$lambda, cvm = cv.lasso$cvm)
    while(lambda_tmp == "NOPE"){
      if(max_f %in% t$nzero){
        lambda_tmp = max(filter(t, nzero==max_f)$lambda)
      }else{
        max_f = max_f-1
      }
    }
    f = coef(cv.lasso, s=lambda_tmp)@Dimnames[[1]][coef(cv.lasso, s = lambda_tmp)@i+1] #Extract features which are kept
  }else{
    lambda_tmp = cv.lasso$lambda.1se
  }
  
#  f = f[2:length(f)]
  f
  
  tmp_coeffs <- coef(cv.lasso, s = lambda_tmp)
  tmp_coeffs = data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
  tmp_coeffs$run = i
  tmp_coeffs$deviance = median(cv.lasso$cvm)
  plot(cv.lasso)
  
  # #Get F1 score from the training data
  pred = as.data.frame(
    glmnet::predict.glmnet(
      object = cv.lasso$glmnet.fit,
      s = lambda_tmp,
      newx = x,
      type="class"
    )
  )
  # 
  colnames(pred) = "Pred"
  pred$"Obs" = y
  pred$Obs[pred$Obs == 1] = "LOF"
  pred$Obs[pred$Obs == 0] = "WT"
  cp_train = cutpointr(pred, x = Pred, class = Obs, pos_class = "LOF",,
                       method = maximize_metric, metric = F1_score,  direction = ">=")
  
  # ggplot(pred, aes(x = Obs, y = Pred))+
  #   geom_boxplot(width = 0.6, size = 0.2, outlier.shape = NA)+
  #   geom_jitter(width= 0.2, size = 0.2)+
  #   geom_hline(yintercept = cp_train$optimal_cutpoint, lty = 2)
  
  #plot_roc(cp_train )
  tmp_coeffs$Gene = paste(gene, collapse=', ' )
  tmp_coeffs$Cohort = cohort
  tmp_coeffs$Study = study
  tmp_coeffs$F1 = cp_train$F1_score
  tmp_coeffs$AUC =   cp_train$AUC
  tmp_coeffs$Sensitivity =   cp_train$sensitivity
  tmp_coeffs$Specificity =   cp_train$specificity
  tmp_coeffs$diff_means = mean(pred$Pred[pred$Obs == "LOF"]) - mean(pred$Pred[pred$Obs == "WT"])

  tmp_coeffs$n_train = paste(table(y)[2], "/", length(y), sep ="")
  
  
  #Test on test data if needed
  if(test == T){
    
   # xnew = dplyr::select(test_data, -(cols_to_remove))
   # xnew = as.matrix(xnew)
    ynew = test_data$y  
    
  pred = as.data.frame(
    glmnet::predict.glmnet(
      object = cv.lasso$glmnet.fit,
      s = lambda_tmp,
     # newx = xnew,
      newx = as.matrix(dplyr::select(test_data, -c( "Donor_ID", "Sample_ID", "y" ,"Study", "primaryTumorLocation", "Mutation_level"))),
      type="class"
    )
  )
  # 
  colnames(pred) = "Pred"
  pred$"Obs" = ynew
  pred$Obs[pred$Obs == 1] = "LOF"
  pred$Obs[pred$Obs == 0] = "WT"
  cp_train = cutpointr(pred, x = Pred, class = Obs, pos_class = "LOF",,
                       method = maximize_metric, metric = F1_score,  direction = ">=")
  
  tmp_coeffs$F1_test = cp_train$F1_score
  tmp_coeffs$AUC_test =   cp_train$AUC
  tmp_coeffs$Sensitivity_test =   cp_train$sensitivity
  tmp_coeffs$Specificity_test =   cp_train$specificity
  tmp_coeffs$diff_means_test = mean(pred$Pred[pred$Obs == "LOF"]) - mean(pred$Pred[pred$Obs == "WT"])
  tmp_coeffs$n_test = paste(table(ynew)[2], "/", length(ynew), sep ="")
  }else{
    tmp_coeffs$F1_test = NA
    tmp_coeffs$AUC_test =   NA
    tmp_coeffs$Sensitivity_test =   NA
    tmp_coeffs$Specificity_test =  NA
    tmp_coeffs$n_test = NA
  }
    
  
  if(i == 1){
    out = tmp_coeffs
  }else{
    out = bind_rows(out, tmp_coeffs)
  }
  print(i)
  
  
  
  #Plot the predictions
  # #cp_train$optimal_cutpoint
  # cp_train$F1_score
  # #plot_precision_recall(cp_train, display_cutpoint = T)
  # plot_sensitivity_specificity(cp_train)
  # plot_roc(cp_train)
  # 
  # ggplot(pred, aes(x = as.factor(y), y = Pred))+
  #   geom_boxplot()+
  #   geom_jitter(width = 0.2, size = 0.2)+
  #   geom_hline(yintercept = cp_train$optimal_cutpoint)
  
  }
  out = out%>%
    group_by(name)%>%
    mutate(n = n())%>%
    ungroup()
  return(out)
}

# BRCA = lassoLoop(train_data = train_data, nIterations = 1, nfold = 5, Null = F, 
#                  test_data = test_data, cols_to_remove = cols_to_remove, test = T,
#                  gene = gene, cohort = cohort, study = trainStudy) #One iteration is plenty, it doesnt change


# 
# tt = lassoLoop(train_data = train_data, nIterations = 100, nfold = 5, Null = T)
# #write.table(tt, paste("Feedback_from_reviewers/Split_by_dataset/NUll_models/",gene,"_",cohort,"_",trainStudy,".tsv", sep =""), col.names = T, row.names = F, quote = F)
# 
# t2 = tt%>%
#   distinct(run, .keep_all = T)
# #hist(t2$AUC, breaks = 10)
# ggplot(t2, aes(x = deviance))+
#   geom_histogram(bins = 10)+
#   geom_vline(xintercept = BRCA$deviance)
# 
# median(tt$AUC)
# hist(tt$F1)
# BRCA$F1
# 
# t2 = tt%>%
#   group_by(name)%>%
#   mutate(n = n())%>%
#   distinct(name, n)





