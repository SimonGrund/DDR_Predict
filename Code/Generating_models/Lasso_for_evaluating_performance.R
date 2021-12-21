library(PRROC)

lassoLoopNestedCv = function(train_data = train_data, cols_to_remove = cols_to_remove, gene = gene, cohort = cohort, study = trainStudy){
  
  #train_data = train_data; test_data = test_data; nIterations = 10; nfold = 5; Null = F; test = F;cols_to_remove = cols_to_remove; gene = gene; cohort = cohort; study = trainStudy
  train_data$ID = 1:nrow(train_data)
  y = train_data$y
  #Create subsets of data for validation
  nfolds = sum(y==1)
  if(nfolds > 30){nfolds = 30}
  nfold <- nfolds
    #sum(y == 1)
  # assign folds evenly using the modulus operator
  fold0 <- sample.int(sum(y==0)) %% nfold
  fold1 <- sample.int(sum(y==1)) %% nfold
  foldid <- numeric(length(y))
  foldid[y==0] <- fold0
  foldid[y==1] <- fold1
  foldid <- foldid + 1
  train_data$fold = foldid
  
  for(i in 1:nfolds){
    
    Validation = dplyr::filter(train_data, fold == i) 
    Training = dplyr::filter(train_data, fold != i) 
    x = dplyr::select(Training, -c(cols_to_remove, "fold", "ID"))
    x = as.matrix(x)
    y = Training$y
    Validation_x = dplyr::select(Validation, -c(cols_to_remove, "fold", "ID"))
    Validation_x = as.matrix(Validation_x)
    Validation_y = Validation$y
    
    weights = Training$weights  
      
  #Start the Lasso
  set.seed(i)
  nfold <- 5
  # assign folds evenly using the modulus operator
  fold0 <- sample.int(sum(y==0)) %% nfold
  fold1 <- sample.int(sum(y==1)) %% nfold
  foldid <- numeric(length(y))
  foldid[y==0] <- fold0
  foldid[y==1] <- fold1
  foldid <- foldid + 1
  
  #  tryCatch({
  cv.lasso <- glmnet::cv.glmnet(
    x = x, 
    y = factor(y), 
    alpha = 1, 
    family = "binomial", 
    type.measure="deviance",   
    intercept = T,
    foldid = foldid,
    #nfolds = nrow(x), #Kfold cross validation, performs worse!
    weights = Training$weights, 
    #lower.limits= 0,
    standardize = F) #Use CV Lasso to estimate best model and return it
  # }, warning=function(cond){e_message <<- "Lasso doesn't converge";next})
  #plot(cv.lasso)
  #We use lambda = 1SE
 # summary(cv.lasso$cvm) #See AUC scores of the runs
  
  # Predict on the validation data
  # pred = as.data.frame( #09-12-2021 , If i could go back I would use this one instead. The glmnet version probs seems buggy, although quite comparable results
  #   predict(
  #     object = cv.lasso$glmnet.fit,
  #     s = cv.lasso$lambda.1se,
  #     newx = Validation_x,
  #     type="response"
  #   )
  # )
  pred = as.data.frame(
    glmnet::predict.glmnet(
      object = cv.lasso$glmnet.fit,
      s = cv.lasso$lambda.1se,
      newx = Validation_x,
      type="response"
    )
  )
  hist(pred$s1)
  
  colnames(pred) = "Posterior_probability"
  pred$ID = Validation$ID
  pred$Obs = Validation$y
  pred$Sample_ID = Validation$Sample_ID
  
  if(i == 1){
    out_pred = pred
  }else{
    out_pred = bind_rows(out_pred, pred)
  }
  
  print(i)
  }

  #Evaluate performance
  fg <- out_pred$Posterior_probability[out_pred$Obs == 1]
  bg <- out_pred$Posterior_probability[out_pred$Obs == 0]
  PRAUC(out_pred$Posterior_probability, y_true = out_pred$Obs)
  pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T) #If fg and bg turned around it is much stronger! Because you add the large rate of the background class :D!
  #plot(pr)
  
  
  pr_AUC = pr$auc.integral
  baseline = sum(out_pred$Obs == 1)/nrow(out_pred)

  # ggplot(train_data, aes(x = as.factor(y), y = Posterior_probability))+
  #   geom_boxplot(width = 0.6, size = 0.2, outlier.shape = NA)+
  #   geom_jitter(width= 0.2, size = 0.2)+
  #   geom_hline(yintercept = cp_train$optimal_cutpoint, lty = 2)

  #Export performance metrics
  y = train_data$y
  out = data.frame(
    Gene = gene,
    Cohort = cohort,
    Study = study,
    AUC = roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)$auc,
    pr_AUC = pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)$auc.integral,
    baseline = baseline,
    n_train = paste(table(y)[2], "/", length(y), sep ="")
    )

  return(out)
}

#set.seed(1)
#BRCA = lassoLoopNestedCv(train_data = train_data, cols_to_remove = cols_to_remove,
   #              gene = gene, cohort = cohort, study = trainStudy) #One iteration is plenty, it doesnt change
#BRCA








