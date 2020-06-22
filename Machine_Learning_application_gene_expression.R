####################################################################
### Aina Montalban - SCRIPT
### Machine Learning model based on gene expression signatures.

### Description:
### Predictor that is able to classify samples responders and non-responders
### Using 3-repeats 10-fold cross-validation
### Uses Random Forest as classification algorithm
### Random forest function should be loaded from the script "customRF.R", 
### an extensive of caret which incorporates hyperparamenter tuning. 

### For the script you should have 4 files:
### a) samplesheet
### b) gene expression dataset
##################################


# Load packages
library(caret)
library(minfi)
library(limma)
library(randomForest)
library(plotROC)
library(reshape2)
library(gridExtra)


source("customRF.R") # we load the Random Forest function that is able to do hyperparamenter tuning
source("/mnt/ElRaid/amontalban/PROJECTS/2020_AMONTALBAN_PDX_Livio/code/FS_method.R") # feature selection methods
# from this script we import the run_fs_methods function, where we can specify which feature selection method we want


# Directories
stem_path <- "/mnt/ElRaid/amontalban/PROJECTS/2020_AMONTALBAN_PDX_Livio" 
ML_DIR <- file.path(stem_path, "aina/results/ML")
ML_RESULTS <- file.path(stem_path, "aina/results/ML")
RESULTS_DIR <- file.path(stem_path, "aina/results")
setwd(RESULTS_DIR)

DATA_DIR <- file.path(
  stem_path, 
  "data")


# Directories
stem_path <- "/mnt/ElRaid/amontalban/PROJECTS/2020_AMONTALBAN_PDX_Livio" 
RESULTS_DIR <- file.path(stem_path, "aina/results")
setwd(RESULTS_DIR)
ML_RESULTS <- file.path(stem_path, "aina/results/mRNA_results")

# Load targets
#"data")
targets <- read.csv(file.path(DATA_DIR, "sample_sheet_PDX.csv"), skip = 7)
dim(targets)
rownames(targets) <- targets$genealogy..LMX..Xenograft.from.liver.mets.
dim(targets)
same_ind <- c( "CRC0356LMX","CRC0307LMX" ,"CRC0161LMX" ,"CRC0358LMX" ,"CRC0262LMX", "CRC0477LMX" ,"CRC0081LMX", "CRC0029LMX", "CRC0080LMX",
               "CRC0144LMX", "CRC0306LMX")
targets <- targets[targets$Sample_Group != "",]
rownames(targets) <- targets$genealogy..LMX..Xenograft.from.liver.mets.
targets <- targets[!rownames(targets) %in% same_ind, ]
dim(targets)

gst <- getGEO("GSE76402", AnnotGPL = TRUE)
gset <- gst[[1]]
ss <- pData(gst[[1]])


subset_target <- targets[,c("Sample_Name", "Sample_Group", "cgenalogy")] #retrieve sample name and group from methylation samplesheet
subset_target$ID <- gsub("LMX","" ,rownames(subset_target)) # remove LMX to match sample names
# subset_target$ID <- paste(subset_target$ID, 
#                           substr(x=subset_target$cgenalogy, start = 12, stop = 12), 
#                           sep = "")

ss$ID <- ss$`case_unique_id:ch1` # create column ID in gene expression samplesheet
#ss$ID <- paste(ss$`case_unique_id:ch1`, ss$`sampling region:ch1`, sep = "") # create column ID in gene expression samplesheet
ss.merged <- merge(subset_target, ss, by="ID") # merge the two samplesheets by ID variable


dim(ss.merged)
length(ss.merged$geo_accession)

# Expression set
eset <- log2(exprs(gset))
eset.merged <- eset[, ss.merged$geo_accession]
colnames(eset.merged) <- ss.merged$ID
rownames(ss.merged) <- ss.merged$ID
dim(eset.merged)


unique_id <- unique(ss.merged$ID) 
ss.unique <- ss.merged[!duplicated(ss.merged$ID),] # select unique samplenames
eset_avg <- t(rowsum(t(eset.merged), colnames(eset.merged)) / c(table(colnames(eset.merged)))) # replicate averaging
eset_avg <- as.data.frame(t(eset_avg))
eset_avg$response <- droplevels(as.factor(ss.unique$Sample_Group))

save(eset_avg, file = file.path(RESULTS_DIR, "eset_avg.RData"))

# Define variables
#== 10-fold cross validation with 3-repeats
nfolds <- 10
nrepeats <- 3
#== Feature selection method: hybrid algorithm using BORUTA and limma
fs_method <- "BORUTA" 

#== Vectors to save the classification metrics in each fold
acc_train_rf <- c()
acc_test_rf <- c()
kappa_test_rf <- c()
sensitivity_test_rf <- c()
specificity_test_rf <- c()

#== List to select the best hyperparameters and the best features selected
best_hyperparameters <- list()
features_selected <- list()

  
# create 80%/20% for training and validation datasets
set.seed(3)
train_index <- createDataPartition(eset_avg$response, p=0.80, list=FALSE)
validation <- eset_avg[-train_index,]
training <- eset_avg[train_index,]
y_training <- droplevels(as.factor(ss.unique$Sample_Group[train_index]))
y_validation <- droplevels(ss.unique$Sample_Group[-train_index])
  
table(y_training) # to know the number of responders and non responders in the training set
table(y_test) # to know the number of responders and non responders in the test set
training$response <- y_training
validation$response <- y_validation
  
index <- 1 #variable to help saving the classification performances
for (j in 1:nrepeats){
    print(paste("Repeat", j))
    set.seed(123)
    folds <- createFolds(training$response, k=nfolds)
    # Inside each fold we include:
    # Feature Selection 
    # Hyperparameters optimization
    for (i in 1:nfolds){
      
      print(paste("Fold", i))
      print(paste("Splitting the data..."))
      
      # Split data into training and test sets
      test <- training[folds[[i]],]
      train <- training[unlist(folds[-i]),]
      # save response variables
      y_test <- test$response
      y_train <- train$response
      
      print(paste("Running feature selection..."))
      # In each fold, perform feature selection
      features  <- run_fs_methods(method=fs_method, num_features=20, trainMatrix=train, y_train)
      # Save features selected to make the model 
      train_fs <- train[, features]
      train_fs$response <- y_train
      # Save selected features for test
      test_fs <- test[, features]
      
      # Run classification algorithm: random Forest
      print(paste("Evaluating RF.."))
      tunegrid <- expand.grid(.mtry=c(2:15),.ntree=c(1000,1500))
      set.seed(123)
      fit.rf <- train(response~., data=train_fs, method=customRF,  metric="Accuracy", tuneGrid=tunegrid) # run RF
      pred_rf_test <- predict(fit.rf, newdata = test_fs, type = "raw") # predict on test data
      acc_test_rf[[index]] <- confusionMatrix(as.factor(pred_rf_test), y_test)$overall["Accuracy"] # calculate accuracy
      kappa_test_rf[[index]] <- confusionMatrix(pred_rf_test, y_test)$overall["Kappa"] # calculate kappa
      sensitivity_test_rf[[index]] <- confusionMatrix(pred_rf_test, y_test)$byClass["Sensitivity"] # calculate sensitivity
      specificity_test_rf[[index]] <- confusionMatrix(pred_rf_test, y_test)$byClass["Specificity"] # calculate specificity
      print(paste("Completed RF algorithm"))
      
      features_selected[[index]] <- features # keep features selected for each splot
      
      best_hyperparameters[[index]] <- fit.rf$bestTune # keep hyperparameeters  for each split 
      
      index <- index + 1
    }
  }
  
write.csv(acc_test_rf, file = file.path(RESULTS_DIR, "acc_folds_exp.csv"))
write.csv(kappa_test_rf, file = file.path(RESULTS_DIR, "kappa_folds_exp.csv"))
write.csv(sensitivity_test_rf, file = file.path(RESULTS_DIR, "sensitivity_folds_exp.csv"))
write.csv(specificity_test_rf, file = file.path(RESULTS_DIR, "specificity_folds_exp.csv"))
  
  
# Prediction on a validation set
hyper <- unname(best_hyperparameters[[which.max(acc_test_rf)]])
best_features <- unlist(features_selected[which.max(acc_test_rf)])
  
fs_train <- training[, best_features]
fs_train$response <- y_training
  
  
# Create the final standalone model using all training data
set.seed(7)
final_model <- randomForest(response~., data=fs_train, mtry=unlist(hyper[1]), ntree=unlist(hyper[2]))
save(final_model, file = file.path(RESULTS_DIR, "final_model_exp.rds")) # save to model for later use
  
  
# Make a predictions on "new data" using the final model
final_predictions <- predict(final_model, validation)
acc <- unname(confusionMatrix(unname(final_predictions), unname(y_validation))$overall["Accuracy"])
kappa <- unname(confusionMatrix(unname(final_predictions), unname(y_validation))$overall["Kappa"])
sensitivity <- unname(confusionMatrix(unname(final_predictions), unname(y_validation))$byClass["Sensitivity"])
specificity <- unname(confusionMatrix(unname(final_predictions), unname(y_validation))$byClass["Specificity"])
  
metrics <- cbind(acc, kappa, sensitivity, specificity) # join accuracy, kappa, sensitivity and specificity
write.csv(metrics, file = file.path(RESULTS_DIR, "metrics.csv")) # save metric on the new data

metrics <- cbind(acc, kappa, sensitivity, specificity)
write.csv(metrics, file = file.path(RESULTS_DIR, "metrics_exp.csv"))

save(final_model, file = file.path(RESULTS_DIR, "final_model_exp.rds"))

# # Create dataframe for ROC curve
ROC_df <- data.frame(row.names = rownames(validation))
ROC_df$D <- unclass(y_validation) - 1
ROC_df$D.str <- y_validation

rf.prob <- predict(final_model, newdata=validation, type="prob")
ROC_df$M_RF  <- rf.prob[,2]
rocRF <- ggplot(ROC_df, aes(d = D, m = M_RF)) + geom_roc(n.cuts = 0) + style_roc()
rocRF <- rocRF + ggtitle("ROC curve") +
  annotate("text", x = .75, y = .25,
           label = paste("AUC =", round(calc_auc(rocRF)$AUC, 3)))
ggsave(file.path(RESULTS_DIR, "Final_ROC_curve_exp.png"))

# Annotate CpGs
annot <- read.delim(file = file.path(RESULTS_DIR, "GPL10558_HumanHT-12_V4_0_R1_15002873_B.txt"), skip = 8)
annot2 <- read.delim(file = file.path(RESULTS_DIR, "GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt"), skip = 8)

imp_var <- varImp(final_model)
imp_var$Probe_Id <- rownames(imp_var)
imp_var_ordered <- imp_var[order(varImp(final_model), decreasing = T),]

probeList <- imp_var_ordered$Probe_Id
annSub <- annot2[annot2$Probe_Id %in% probeList, c("Probe_Id", "Symbol", "Definition")]
imp_var_sym <- merge(imp_var_ordered, annSub, by="Probe_Id")
imp_var_final <- imp_var_sym[order(imp_var_sym$Overall, decreasing = T),]
imp_var_final$Symbol
write.csv(imp_var_final, file = file.path(RESULTS_DIR, "imp_pred_final_exp.csv"))
