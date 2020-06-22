####################################################################
### Aina Montalban - SCRIPT
### Feature selection methods function

### Description:
### Feature selection methods used to select the optimal subset of featuers
### Filter methods: LIMMA and ANOVA
### Hybrid method: w/ filter methods and BORUTA / L1norm / RFS
### Regions differentially methylated


### For the script you should have 4 files:
### a) training matrix or data frame with the response variable
### b) number of features to select
### c) vector with the response variable
##################################

# Load packages
library(minfi)
library(limma)
library(caret)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(doMC)
require(parallel)
require(stats)
annot <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)


# LIMMA feature selection: it performs modeterated t-test in the training matrix

FS_limma <- function(trainMatrix, num_features = 1000, contr, y_train){
  my_factor <- factor(trainMatrix$response)
  design <- model.matrix(~0+my_factor, data = trainMatrix$response)
  colnames(design) <- levels(my_factor)
  
  fit <- lmFit(t(trainMatrix[, -c(ncol(trainMatrix))]), design)
  contMatrix <- makeContrasts(contrasts = contr, levels=design)
  
  fit2 <- contrasts.fit(fit, contMatrix)
  fit2 <- eBayes(fit2)
  
  print(summary(decideTests(fit2)))
  top <- topTable(fit2, coef = 1, number = Inf, sort.by="P")
  cpgs <- rownames(top)[1:num_features]
  return(cpgs)
}


# ANOVA feature selection
FS_anova <- function(trainMatrix, num_features, y_train){
  y <- train$response
  anova_test <- mclapply(train[, -ncol(train)], function(x){
    cpg_data <- data.frame(x, y)
    print(head(x))
    print(head(y))
    
    test <- aov(x ~ y, data = cpg_data)
    results <- summary(test)
    print(results)
    fval <- results[[1]]$`F value`[1]
    pval <- results[[1]]$`Pr(>F)`[1]
    return(pval)
  }, mc.cores = 4)
  
  cpgs_anova <- sort(unlist(anova_test))
  print(head(cpgs_anova))
  cpgs_anova <- names(cpgs_anova)[1:num_features]
  
}

# BORUTA feature selection
FS_boruta <- function(training, num_features, y_response){
  boruta.train <- Boruta::Boruta(x=training, y=y_response, doTrace=0)
  boruta.res <- Boruta::attStats(boruta.train)
  
  if ("Tentative" %in% boruta.res$decision) {
    final.boruta <- Boruta::TentativeRoughFix(boruta.train)
    boruta.res <- Boruta::attStats(final.boruta)
  }
  # Sort features by mean Importance
  boruta.res <- boruta.res[order(boruta.res$meanImp,
                                 decreasing = TRUE), ]
  # Select all confirmed important attributes, ordered by meanImp.
  cpgs_boruta <- rownames(boruta.res)[boruta.res$decision=="Confirmed"]
  if (length(cpgs_boruta) <= num_features){
    return(cpgs_boruta)
  } else{
  return(cpgs_boruta[1:num_features])
  }
}

# Differentially methylated regions feature selection function
FS_dmr <- function(training, num_regions){
  my_factor <- factor(training$response)
  design <- model.matrix(~0+my_factor, data = training$response)
  colnames(design) <- levels(my_factor)
  
  contMatrix <- makeContrasts(contrasts="NonResponder-Responder", levels=design)
  
  myAnnotation <- DMRcate::cpg.annotate(object=t(training[,-ncol(training)]), 
                                       datatype="array", what="M", analysis.type="differential", 
                                       design=design, contrasts=TRUE, cont.matrix=contMatrix, 
                                       coef="NonResponder-Responder", arraytype="EPIC", fdr=0.05)
  DMR <- DMRcate::dmrcate(myAnnotation, lambda=1000, C=2)
  
  # convert regions to annotated genomic ranges
  results.ranges.dmr <- DMRcate::extractRanges(DMR, genome="hg19")
  
  gen <- "hg19"
  dmr_idxs <- 1:num_regions
  
  annot <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  annSub <- annot[match(colnames(training[,-ncol(training)]),annot$Name), c(1:4,12:19, 22:24, 31:ncol(annot))]
  
  
  annSubOrd <- annSub[order(annSub$chr, annSub$pos),]
  traning_Ord <- training[, match(annSubOrd$Name, colnames(training))]
  
  dataset_DMR <- data.frame(matrix(ncol = 0, nrow = nrow(training)))
  region_cpg <- list()
  # Create DMR regions dataset
  for (dmrIndex in dmr_idxs){
    print(paste("Computing region ", dmrIndex, "..." , sep=""))
    # extract chromosome number and location from DMR results 
    chrom <- as.character(seqnames(results.ranges.dmr[dmrIndex]))
    start <- as.numeric(start(results.ranges.dmr[dmrIndex]))
    end <- as.numeric(end(results.ranges.dmr[dmrIndex]))
    
    # create genomic ranges object from methylation data
    cpgData <- GRanges(seqnames=Rle(annSubOrd$chr),
                       ranges=IRanges(start=annSubOrd$pos, end=annSubOrd$pos),
                       strand=Rle(rep("*",nrow(annSubOrd))),
                       betas=t(traning_Ord))
    
    # extract data on CpGs in DMR
    cpgData <- subsetByOverlaps(cpgData, results.ranges.dmr[dmrIndex])
    df_cpgData <- as.data.frame(cpgData)
    cpg <- annSubOrd$Name[match(df_cpgData$start, annSubOrd$pos)] 
    df_cpgData <- df_cpgData[, 6:ncol(df_cpgData)]
    dataset_DMR$name_region <- unname(colMeans(df_cpgData))
    colnames(dataset_DMR)[dmrIndex] <- paste("region", dmrIndex, sep = ".")
    
    region_cpg[[dmrIndex]] <- cpg
    names(region_cpg)[dmrIndex] <- paste("region", dmrIndex, sep = ".")
  }
  return(list(regions2cpg=region_cpg, data=dataset_DMR))
  
}


# L1-norm feature selection method
FS_L1norm <- function(training, num_features, y_response){
  set.seed(3)
  training$response <- y_response
  x <- model.matrix(as.factor(response)~., training)[, -1]
  #set.seed(3)
  registerDoMC(cores = 4)
  cv.lasso <- cv.glmnet(x, training$response, family='binomial', alpha=1, parallel=TRUE, type.measure='auc')
  plot(cv.lasso)
  
  coeff <- as.matrix(coef(cv.lasso, s=cv.lasso$lambda.min))
  res <- coeff[coeff[, 1] != 0, ]
  res <- res[order(res)]
  cpg_L1norm <- names(res)
  cpg_L1norm <- cpg_L1norm[-which(cpg_L1norm == "(Intercept)")]
  if (length(cpg_L1norm) <= num_features){
    return(cpg_L1norm)
  } else{
    return(cpg_L1norm[1:num_features])
  }
}


# Recursive Feature Selection method 
FS_RFS <- function(training, num_features, y_response){
        ctrl <- rfeControl(method = "cv", 
                   verbose = TRUE, functions = rfFuncs, allowParallel = TRUE)
        registerDoMC(cores = 4)
        rfRFE <- rfe(x=training, y=y_response, 
             sizes = c(30:100, 200, 250, 300),
             metric="Accuracy", rfeControl = ctrl, ntree=1000)
        return(predictors(rfRFE)[1:num_features])
}


## Global function to choose the feature selection method:
run_fs_methods <- function(method, num_features, trainMatrix, y_train){
  if (method == "LIMMA"){
    cpgs <- FS_limma(trainMatrix, num_features , "NonResponder-Responder", y_train)
    return(cpgs)
  } else if (method == "ANOVA"){
    cpgs <- FS_anova(trainMatrix, num_features, y_train)
    return(cpgs)
  } else if (method == "BORUTA"){
    diff_cpgs <- FS_limma(trainMatrix, num_features = 1000 , "NonResponder-Responder", y_train)
    cpgs <- FS_boruta(trainMatrix[, diff_cpgs], num_features, y_train)
  } else if (method == "DMR"){
    cpgs <- FS_dmr(trainMatrix, num_features)
    return(cpgs)
  } else if (method == "L1norm"){
    diff_cpgs <- FS_limma(trainMatrix, num_features = 1000 , "NonResponder-Responder", y_train)
    cpgs <- FS_L1norm(trainMatrix[, diff_cpgs], num_features, y_train)
    return(cpgs)
  } else if (method == "RFS"){
    diff_cpgs <- FS_limma(trainMatrix, num_features = 200 , "NonResponder-Responder", y_train)
    cpgs <- FS_RFS(trainMatrix[, diff_cpgs], num_features, y_train)
    return(cpgs)
  }
}



