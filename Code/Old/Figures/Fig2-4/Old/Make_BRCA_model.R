library(data.table)
library(tidyverse)
library(caret)
setwd("/Volumes/GenomeDK/PCAWG/simong_SV/DDR_project_Jan_2021//")

d_lof = fread("Data/LOF_all.tsv")
all_samples = fread("Data/all_donors.tsv")
features = fread("Data/COMBINED_featureCounts.tsv")
colnames(features) = str_replace_all(colnames(features), "-",".")
colnames(features) = str_replace_all(colnames(features), ">","_._")
d_lof = filter(d_lof, Sample_ID %in% features$Sample_ID)
all_samples = filter(all_samples, Sample_ID %in% features$Sample_ID)

#Filter all data
d = filter(d_lof, GENE %in% c("BRCA2", "BRCA1"), primaryTumorLocation %in% c("Breast", "Prostate", "Ovary", "Pancreas"))
d2 = anti_join(all_samples, d)%>%
  mutate(LOF = 0, GENE = "BRCA", Germline = F, Double = F)
d = bind_rows(d, d2)%>%
  filter( primaryTumorLocation %in% c("Breast", "Prostate", "Ovary", "Pancreas"))%>%
  dplyr::select(-cancertype)

d = left_join(d, features)

#Get a one-value per row
table(is.na(d))

d = d%>%
  group_by(Sample_ID)%>%
  mutate(LOF = max(LOF),
         Double = ifelse(sum(Double == T)>0,T,F),
         Germline = ifelse(sum(Germline == T)>0,T,F),
  )%>%
  ungroup()%>%
  distinct(Sample_ID, .keep_all = T)
d = mutate(d, y = ifelse(LOF >= 6, 1,0))
ds = d
d = filter(d, !LOF %in% c(2,3,4,5))

#Run the lasso
train_data = d

w = table(train_data$y)
w1 = w[1]/nrow(train_data)
w0 = w[2]/nrow(train_data)
train_data = mutate(train_data,
                    y = as.factor(y),
                    weights = ifelse(y == 1, w1, w0),
                    weights = ifelse(Germline == T, (weights * 1.2), weights),
                    weights = ifelse(Double == T, (weights * 1.2), weights),
                    weights = ifelse(LOF == 7, (weights * 1.2), weights)
)%>% #Add 20% weight if germline,another 20 if clinvar
  dplyr::select(Donor_ID, Sample_ID, LOF, GENE, y, Double, weights, Germline, LOF, Study, primaryTumorLocation, everything())

cols_to_remove = c("Donor_ID", "Sample_ID", "LOF","Double", "GENE", "y", "weights", "Germline", "LOF", "Study", "primaryTumorLocation")

####
# CV.LASSO
####
x = as.matrix(dplyr::select(train_data, -(cols_to_remove)))
#x <- sapply( x, as.numeric )
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
                      intercept = F,
                      foldid = foldid,
                      nfolds = nrow(x),
                      weights = train_data$weights,
                      standardize = T) #Use CV Lasso to estimate best model and return it
# }, warning=function(cond){e_message <<- "Lasso doesn't converge";next})

plot(cv.lasso)
max_f = 8
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

f


#formula  = as.formula(paste("y~primaryTumorLocation*(",paste(f, collapse = "+"),")", sep = ""))
formula  = as.formula(paste("y~primaryTumorLocation*(",paste(f, collapse = "+"),")", sep = ""))

#Make GLM
train_data2 = as.data.frame(train_data)%>%
  mutate(
    y = ifelse(y == 1, "LOF", "WT")
  )

tmp = train_data2
tmp[,f] = as.numeric(scale(log(tmp[,f]+0.1))) #Scale it allll together :D
tmp$y= factor(tmp$y, levels = c("WT", "LOF")) #Set the WT as the 'positive' class to get estimates for the transition into WT
train_data2 = tmp

train_data2$y = factor(train_data2$y, levels = c("WT", "LOF"))
table(train_data2$y)
c = caret::train(formula, 
                 method = "glm", 
                 data = train_data2, 
                 trControl = trainControl(method = "LOOCV", classProbs = TRUE),
                 weights = train_data2$weights,
                 family = "quasibinomial"
)
summary(c)
hist(c$pred$LOF)
#plot(c)
summary(c)
#as.matrix(train_data2[,f])
#p = predict(c, newdata = as.matrix(train_data2[,f]))
#hist(p$LOF)

t = as.data.frame(summary(c)$coefficients)
Features = data.frame(Gene = "BRCA1/2", 
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

p2 = ggplot(Features, aes(x= Estimate, y = Feature))+
  geom_errorbar(aes(xmin = Estimate - stderror, xmax = Estimate + stderror))+
  facet_wrap(~paste(Gene, sep = ", "), scales = "free", nrow = 6)+
  theme_bw(base_size = 8, base_family = "Arial") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.9, 0.1),
        panel.border = element_rect(color = "black", fill = NA),
        axis.line = element_line(color = NA)) +
  geom_vline(xintercept = 0, lty = 2, col = "grey")

p2

hist(c$pred$LOF)
train_set = c$pred
train_set$mutation_level = train_data2$LOF 
train_set$"Source" = "Training"
train_set = dplyr::select(train_set, -weights, -rowIndex, -parameter)
train_set$Sample_ID = train_data2$Sample_ID
train_set$CT = train_data$primaryTumorLocation

cp_train = cutpointr(train_set, x = LOF, class = obs, pos_class = "LOF",
          method = maximize_metric, metric = sum_sens_spec,  direction = ">=", subgroup = CT)

plot_roc(cp_train)
cp_train$AUC 

saveRDS(object = c, file = paste("Figures/Figure_3/BRCA_model.Rdata" ,sep = ""))

#Run model on all data—also VUS

ds2 = as.data.frame(ds)%>%
  mutate(
    y = ifelse(y == 1, "LOF", "WT")
  )

tmp = ds2
tmp[,f] = as.numeric(scale(log(tmp[,f]+0.1))) #Scale it allll together :D
tmp$y= factor(tmp$y, levels = c("WT", "LOF")) #Set the WT as the 'positive' class to get estimates for the transition into WT
ds2 = tmp

ds2$y = factor(ds2$y, levels = c("WT", "LOF"))
p = caret::predict.train(c, newdata = ds2, type = "prob")
#p = as.data.frame(predict(object = c, newdata = ds2))
#colnames(p) = "LOF"

#Add true status
p$"Sample_ID" = ds$Sample_ID
p$"primaryTumorLocation" = ds$primaryTumorLocation
p$obs = ds2$y
p$Mutation_level = ds2$LOF
#hist(p$LOF)
cutpointr(p, x = LOF, class = obs, pos_class = "LOF",
          method = maximize_metric, metric = sum_sens_spec, direction = ">=", subgroup = primaryTumorLocation)

write.table(x = p, file =  "Figures/Figure_3/BRCA_predicted_data.tsv", sep ="\t", col.names = T, row.names = F)

