####################################################################
### Aina Montalban
### GO terms plots
### Description:
### Intersection between differentially expressed genes (DEGs) and differentially methylated genes (DMPs)
### Find DMP-DEG pairs correlation 
### For the script you should have 4 files:
### a) GO_Biological_Process.csv
### b) GO_Cellular_Component.csv
### c) GO_Molecular_Function.csv
##################################


# Machine Learning model
library(caret)
library(minfi)
library(limma)
library(randomForest)
library(e1071)
library(class)
library(MASS)
library(glmnet)
library(gbm)
library(plotROC)
library(reshape2)
library(gridExtra)
# Directories
stem_path <- "/mnt/ElRaid/amontalban/PROJECTS/2020_AMONTALBAN_PDX_Livio" 
RESULTS_DIR <- file.path(stem_path, "results")
setwd(RESULTS_DIR)
ML_RESULTS <- file.path(stem_path, "aina/results/ML_BVAL_RESULTS")

DATA_DIR <- file.path(
  stem_path, 
  "data")
# Load targets
targets <- read.csv(file.path(DATA_DIR, "sample_sheet_PDX.csv"), skip = 7)

same_ind <- c( "CRC0356LMX","CRC0307LMX" ,"CRC0161LMX" ,"CRC0358LMX" ,"CRC0262LMX", "CRC0477LMX" ,"CRC0081LMX", "CRC0029LMX", "CRC0080LMX",
               "CRC0144LMX", "CRC0306LMX")
targets <- targets[targets$Sample_Group != "",]
rownames(targets) <- targets$genealogy..LMX..Xenograft.from.liver.mets.
targets <- targets[!rownames(targets) %in% same_ind, ]
dim(targets)
# Load b-values and m-values
load("meth.RData")

df_mVals <- as.data.frame(t(mVals))
dim(df_mVals)
df_mVals$response <- droplevels(as.factor(targets$Sample_Group))
# save(df_mVals, file = file.path(RESULTS_DIR, "df_mVals.rds"))


# import feature selection functions
source("/mnt/ElRaid/amontalban/PROJECTS/2020_AMONTALBAN_PDX_Livio/code/FS_method.R")


new_test <- NULL
DMP_2_DMR_test <- function(regionsList, test){
  for (region in 1:length(regionsList)){
    features <- regionsList[[region]]
    new_test <- cbind(new_test, unname(rowMeans(test[, features])))
  }
  colnames(new_test) <- names(regionsList)
  rownames(new_test) <- rownames(test)
  return(new_test)
}



# Variables
nfolds <- 10
acc_train_nb <- c()
acc_train_rf <- c()
acc_train_gbm <- c()
acc_train_elnet <- c()
acc_train_lasso <- c()
acc_train_rda <- c()
acc_train_svm <- c()

acc_test_nb <- c()
acc_test_rf <- c()
acc_test_gbm <-c()
acc_test_elnet <- c()
acc_test_lasso <-c()
acc_test_rda <- c()
acc_test_svm <- c()


kappa_test_nb <- c()
kappa_test_svm <- c()
kappa_test_rf <- c()
kappa_test_gbm <- c()
kappa_test_elnet <- c()
kappa_test_lasso <- c()
kappa_test_rda <- c()


df_cls_test <- data.frame(row.names = c("NB", "RF", "LASSO", "ELNET", "RDA", "GBM", "SVM"))
df_cls_train <- data.frame(row.names = c("NB", "RF", "LASSO", "ELNET", "RDA", "GBM", "SVM"))
kappa_cls_train <- data.frame(row.names = c("NB", "RF", "LASSO", "ELNET", "RDA", "GBM", "SVM"))

fs_method <- "DMR"

print(paste("Using", fs_method))
folds <- createFolds(df_mVals$response, k=nfolds)
for (i in 1:nfolds){
      print(paste("Fold", i))
      print(paste("Splitting the data..."))
      # Split data into training and test sets
      test <- df_bVals[folds[[i]],]
      train <- df_bVals[unlist(folds[-i]),]
      # save response variables
      y_test <- test$response
      y_train <- train$response
      
      
      print(paste("Running feature selection..."))
      # Feature selection
      
      if (fs_method == "DMR"){
        res <- run_fs_methods(method="DMR", num_features=30, trainMatrix=train, y_train)
        
        # Training
        train_fs <- res$data
        train_fs$response <- y_train
        
        # Test
        regList <- res$regions2cpg
        new_test <- NULL
        test_fs <- DMP_2_DMR_test(regList, test)
        
      }else{
        features  <- run_fs_methods(method=fs_method, num_features=20, trainMatrix=train, y_train)
        # Training
        train_fs <- train[, features]
        train_fs$response <- y_train
        
        # Test
        test_fs <- test[, features]
      }

      

      # NB
      print(paste("Evaluating NB.."))
      set.seed(3)
      fit.nb <- train(response~., data = train_fs, method="nb", metric="Accuracy")
      pred_nb_train <- predict(fit.nb, newdata = train_fs)
      acc_train_nb[[i]] <- confusionMatrix(as.factor(pred_nb_train), y_train)$overall["Accuracy"]
      pred_nb_test <- predict(fit.nb, newdata = test_fs, type = "raw")
      acc_test_nb[[i]] <- confusionMatrix(as.factor(pred_nb_test), y_test)$overall["Accuracy"]
      kappa_test_nb[[i]] <- confusionMatrix(as.factor(pred_nb_test), y_test)$overall["Kappa"]
    
      
      
      #  SVM
      print(paste("Evaluating SVM.."))
      set.seed(3)
      fit.svm <- svm(response~., data=train_fs, kernel="linear", scale=FALSE, probability=T)
      pred_svm_train <- predict(fit.svm, newdata = train_fs)
      acc_train_svm[[i]] <- confusionMatrix(as.factor(pred_svm_train), y_train)$overall["Accuracy"]
      predictions.svm <- predict(fit.svm, newdata=test_fs, probability = TRUE)
      acc_test_svm[[i]] <- confusionMatrix(predictions.svm, y_test)$overall["Accuracy"]
      kappa_test_svm[[i]] <- confusionMatrix(predictions.svm, y_test)$overall["Kappa"]
    
      
      # Random Forest
      print(paste("Evaluating RF.."))
      tuneGrid <- expand.grid(.mtry=c(1:15))
      fit.rf <- train(response~., data = train_fs, method="rf", metric="Accuracy", tuneGrid=tuneGrid)
      pred_rf_train <- predict(fit.rf, newdata = train_fs)
      acc_train_rf[[i]] <- confusionMatrix(as.factor(pred_rf_train), y_train)$overall["Accuracy"]
      pred_rf_test <- predict(fit.rf, newdata = test_fs, type = "raw")
      acc_test_rf[[i]] <- confusionMatrix(as.factor(pred_rf_test), y_test)$overall["Accuracy"]
      kappa_test_rf[[i]] <- confusionMatrix(pred_rf_test, y_test)$overall["Kappa"]
      
    
      
      # Regularization methods
      ## Data preparation
      set.seed(3)
      x.train <- model.matrix(as.factor(response)~., train_fs)[, -1]
      test.lasso <- as.data.frame(test_fs)
      test.lasso$response <- y_test
      x.test <- model.matrix(as.factor(response)~., test.lasso)[, -1]
     
      # LASSO
      print(paste("Evaluating LASSSO..."))
      tuneGrid <- expand.grid(alpha = 1,lambda = seq(0.01, 0.03, by = 0.002))
      fit.lasso <- train(response~., data = train_fs, method="glmnet", metric="Accuracy",  tuneGrid = tuneGrid, trControl=caret::trainControl())
      predictions.lasso_train <- predict(fit.lasso, newdata = x.train, type = "raw")
      acc_train_lasso[[i]] <- confusionMatrix(as.factor(predictions.lasso_train), y_train)$overall["Accuracy"]
      predictions.lasso_test <- predict(fit.lasso, newdata = x.test, type = "raw")
      acc_test_lasso[[i]] <- confusionMatrix(as.factor(predictions.lasso_test), y_test)$overall["Accuracy"]
      kappa_test_lasso[[i]] <- confusionMatrix(predictions.lasso_test, y_test)$overall["Kappa"]
      


      # ELNET
      print(paste("Evaluating ELNET..."))
      ## ELASTIC NET ALGORITHM
      tuneGrid <- expand.grid(alpha = seq(0, 1, by=.2),lambda = seq(0.01, 0.03, by = 0.002))
      fit.elnet <- train(response~., data = train_fs, method="glmnet", metric="Accuracy",  tuneGrid = tuneGrid, trControl=caret::trainControl())
      predictions.elnet_train <- predict(fit.elnet, newdata = x.train, type = "raw")
      acc_train_elnet[[i]] <- confusionMatrix(as.factor(predictions.elnet_train), y_train)$overall["Accuracy"]
      predictions.elnet_test <- predict(fit.elnet, newdata =x.test, type = "raw")
      acc_test_elnet[[i]] <- confusionMatrix(as.factor(predictions.elnet_test), y_test)$overall["Accuracy"]
      kappa_test_elnet[[i]] <- confusionMatrix(predictions.elnet_test, y_test)$overall["Kappa"]

      
      # RDA
      print("Running RDA...")
      fit.rda  <- train(response ~ ., data=train_fs, method="rda", metric="Accuracy")
      pred_rda_train <- predict(fit.rda, newdata = train_fs)
      acc_train_rda[[i]] <- confusionMatrix(as.factor(pred_rda_train), y_train)$overall["Accuracy"]
      pred_rda_test <- predict(fit.rda, newdata =test_fs, type = "raw")
      acc_test_rda[[i]] <- confusionMatrix(as.factor(pred_rda_test), y_test)$overall["Accuracy"]
      kappa_test_rda[[i]] <- confusionMatrix(pred_rda_test, y_test)$overall["Kappa"]
      
      
      ### 6 GBM
      print("Running GBM...")
      fit.gbm  <- train(response ~ ., data=train_fs, method="gbm", metric="Accuracy")
      pred_gbm_train <- predict(fit.gbm, newdata = train_fs)
      acc_train_gbm[[i]] <- confusionMatrix(as.factor(pred_gbm_train), y_train)$overall["Accuracy"]
      pred_gbm_test <- predict(fit.gbm, newdata = test_fs, type = "raw")
      acc_test_gbm[[i]] <- confusionMatrix(as.factor(pred_gbm_test), y_test)$overall["Accuracy"]
      kappa_test_gbm[[i]] <- confusionMatrix(predictions.svm, y_test)$overall["Kappa"]
}



df_cls_test <- as.data.frame(rbind(acc_test_nb,  acc_test_rf,
                                   acc_test_elnet, acc_test_lasso, 
                                   acc_test_rda, acc_test_gbm, acc_test_svm))
colnames(df_cls_test) <- c("Fold1", "Fold2", "Fold3" ,"Fold4", "Fold5", 
                           "Fold6", "Fold7", "Fold8" ,"Fold9", "Fold10")
rownames(df_cls_test) <- c("NB", "RF", "ELNET" ,"LASSO", "RDA", "GBM", "SVM")
df_cls_test$method <- rownames(df_cls_test)
write.csv(df_cls_test, file.path(ML_RESULTS, paste("Accuracy_Test",fs_method , ".csv", sep="")))  

kappa_cls_test <- as.data.frame(rbind(kappa_test_nb,  kappa_test_rf,
                                   kappa_test_elnet, kappa_test_lasso, 
                                   kappa_test_rda, kappa_test_gbm, kappa_test_svm))
colnames(kappa_cls_test) <- c("Fold1", "Fold2", "Fold3" ,"Fold4", "Fold5", 
                              "Fold6", "Fold7", "Fold8" ,"Fold9", "Fold10")
rownames(kappa_cls_test) <- c("NB", "RF", "ELNET" ,"LASSO", "RDA", "GBM", "SVM")
kappa_cls_test$method <- rownames(df_cls_test)
write.csv(kappa_cls_test, file.path(ML_RESULTS, paste("Kappa_Test",fs_method , ".csv", sep="")))  


df_cls_train <- as.data.frame(rbind(acc_train_nb,  acc_train_rf,
                                    acc_train_elnet, acc_train_lasso, acc_train_rda, 
                                    acc_train_gbm, acc_train_svm))
colnames(df_cls_train) <- c("Fold1", "Fold2", "Fold3" ,"Fold4", "Fold5", 
                            "Fold6", "Fold7", "Fold8" ,"Fold9", "Fold10")
rownames(df_cls_train) <- c("NB", "RF", "ELNET" ,"LASSO", "RDA", "GBM", "SVM")
write.csv(df_cls_train, file.path(ML_RESULTS, paste("Accuracy_Train",fs_method , ".csv", sep="")))  

      
  

long_df_test <- melt(df_cls_test, "method")
    
ggplot(long_df_test, aes(x=method, y=value, fill=method)) + 
  geom_boxplot(show.legend = FALSE) +
  labs(x="Classifcation Algorithm", y="Accuracy") + 
  scale_fill_brewer(palette = "Dark2") + theme_classic()
ggsave(file.path(ML_RESULTS, paste("Accuracy_Boxplot",fs_method , ".png", sep="")))
        # 


