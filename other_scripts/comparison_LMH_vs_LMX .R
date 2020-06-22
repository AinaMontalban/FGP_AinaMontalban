####################################################################
### Aina Montalban - SCRIPT
### Comparison between LMX-LMH.

### Description:
### DNA methylome comparison between patient-derived xenografts and human liver metastasis

### For the script you should have 4 files:
### a) RAW idat files from both datasets (LMX-LMH)
### b) Samplesheet
##################################
## PACKAGES.
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylationEPICmanifest)
library(RColorBrewer)
library(Rtsne)
library(missMethyl)
library(matrixStats)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
library(gplots)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
#devtools::install_github("jkrijthe/Rtsne")
# install.packages("gplots")
#===============================================================================
### GLOBAL VARIABLES.
stem_path <- "/mnt/ElRaid/amontalban/PROJECTS/2020_AMONTALBAN_PDX_Livio" # dori

# Raw data folder, with the sample sheet and the idat files into subfolders.
IDAT_DIR <- file.path(stem_path, "LMH_vs_LMX_idats")

# Where results will be written. It should exists.
RESULTS_DIR <- file.path(
  stem_path, 
  "aina/results")

# Path where project data (as samplesheet.csv) is stored.
DATA_DIR <- file.path(
  stem_path, 
  "data")


SAMPLESHEET <- "subset_LMX_vs_LMH.csv"

# csv files containing probes to exclude as they are cross-reactive.
POLYMORPHICEPIC <- file.path(
  stem_path,
  "illumina450k_filtering-master/EPIC/13059_2016_1066_MOESM1_ESM.csv")
# Sample sheet column name to use as in differential methylation tests.

class_label <- "Sample_Group"


# Sample sheet column that identify individuals (or any other paired info), for 
# paired tests. If not paired test, let "".
#individual_label <- "Cell_line"
# Contrast to perform (to use in makeContrasts).
contr <- c(
  "HS-Xenograft"
)


autopilot <- FALSE


norm_method <- "Noob"



# Quality filters:
check_sex <- FALSE

genotyping_SNP_QC <- TRUE

# Remove poor quality samples based on mean detection pval.
sample_mean_detP_max <- 0.05

# Remove poor quality samples based on proportion of failed probes.
max_failed_probes <- 0.05

# Decide whether to remove probes failed in any sample or to mark them as NAs
remove_or_NA_probes <- "remove"

# Threshold to consider failed probe.
probe_detP_max <- 0.01

# If your data includes males and females, remove probes on the sex chromosomes
remove_sex_chr <- TRUE

# Remove probes with SNPs at CpG site or at a single nucleotide extension.
remove_SNP_probes <- TRUE

# Exclude cross reactive probes 
remove_cross_reactive <- TRUE

# Perform a variance filtering to filter out bVals under average variance.
perform_var_filtering <- FALSE

# Proportion of most variable bVals to keep, in case variance filtering is TRUE.
# This var can be a proportion in which case this proportion (e.g: 1/3, 0.5, 
# etc) of most variable probes will be selected, or the word "average", in which
# case selected probes will be those >= average probe variance.
var_prop_to_keep <- "average"

#===============================================================================
### FUNCTIONS.
#' Detect Methylation Array platform

detect_array <- function(x) {
  if (x == "IlluminaHumanMethylation450k") {
    return("450k")
  } else if (x == "IlluminaHumanMethylationEPIC") {
    return("EPIC")
  } else {
    return("Unknown")
  }
}

#' Drop CpGs by SNP annotation
#' 
#' This is a custom version of minfi::dropLociWithSnps intented to be used
#' with a beta-values or m-values matrix as input, instead of a 
#' GenomicRatioSet object, despite of it also accepts GenomicRatioSet objects,
#' in which case returns behaves as minfi::dropLociWithSnps and returns a 
#' filtered GenomicRatioSet
#' 
#' @param a_matrix a beta-values or m-values matrix. A GenomicRatioSet object
#'                 is also possible, in which case this function behaves as
#'                 minfi::dropLociWithSnps
#' @param snpAnno a DataFrame object generated with minfi::getSnpInfo function
#'                containing the SNP annotation of the `a_matrix` probes
#' @param snps a character vector with the filtering to be used. 
#'             Default: c("CpG", "SBE")
#' @param maf a numeric value of minor allele frequency to filter out probes.
#'            Default: 0
#' 
#' @details Three filter criteria can be used: "Probe", "CpG" and "SBE" 
#'          correspond the SNPs present inside the probe body, at the CpG 
#'          interrogation and at the single nucleotide extension respectively
#' 
#' @return The same matrix as `a_matrix` but with the SNP probes filtered out.
#'         If a_matrix is a GenomicRatioSet object, this same object class is
#'         returned 
dropLociWithSnps_custom <- function(a_matrix, 
                                    snpAnno, 
                                    snps = c("CpG", "SBE"), 
                                    maf = 0) {
  res <- a_matrix
  to_remove <- character(0)
  # Drop "Probe" SNPs
  if ("Probe" %in% snps) {
    to_remove <- c(to_remove, 
                   rownames(snpAnno)[!(is.na(snpAnno$Probe_rs)) & 
                                       snpAnno$Probe_maf > maf])
  }
  
  # Drop "CpG" SNPs
  if ("CpG" %in% snps) {
    to_remove <- c(to_remove, 
                   rownames(snpAnno)[!(is.na(snpAnno$CpG_rs)) & 
                                       snpAnno$CpG_maf > maf])
  }
  
  # Drop "SBE" SNPs
  if ("SBE" %in% snps) {
    to_remove <- c(to_remove, 
                   rownames(snpAnno)[!(is.na(snpAnno$SBE_rs)) & 
                                       snpAnno$SBE_maf > maf])
  }
  to_keep <- !(rownames(res) %in% to_remove)
  res <- res[to_keep,]
  return(res)
}


#===============================================================================
### MAIN.
################################################################################
## Step 1: data loading
################################################################################
# Read in the sample sheet for the experiment.
targets <- read.metharray.sheet(DATA_DIR, pattern=SAMPLESHEET)

if (!autopilot) {
  print(targets)
}

# Repair Basename col as it is not well done.
targets$Basename <- gsub("_Grn.idat","",list.files(IDAT_DIR, pattern="Grn.idat", full.names = T))

# Enforce Sample_Name as character.
targets$Sample_Name <- as.character(targets$Sample_Name)
print(paste("Full sample sheet rows (samples):", nrow(targets)))


# Read in the raw data from the IDAT files.
rgSet1 <- read.metharray.exp(targets=targets, force = TRUE)
rgSet1


# Detection of Illumina methylation array platform from rgSet.
array_platform <- detect_array(rgSet@annotation[1])
# get the annotation data and probes to exclude.
if (array_platform == "450k") {
  annot <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  probes_to_exclude <- POLYMORPHIC450k
} else if (array_platform == "EPIC") {
  annot <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  probes_to_exclude <- POLYMORPHICEPIC
} else {
  stop("[ERROR] Methylation array platform not detected.")
}
if (!autopilot) {
  print(head(annot))
}
probe_detP_max <- 0.01

################################################################################
## Step 2: QC
################################################################################

snps <- getSnpBeta(rgSet1)
colnames(snps) <- targets$Sample_Name
my_palette <- colorRampPalette(c("darkred", "darkgreen", "darkblue"))(n = 299)
heatmap.2(snps, col = my_palette, trace = "none", margins = c(7, 4.5))

# calculate the detection p-values
detP <- detectionP(rgSet1)

pal <- brewer.pal(8,"Dark2")
mean_detP <- colMeans(detP)
names(mean_detP) <- targets$Sample_Name
barplot(mean_detP, col=pal[factor(targets[, class_label])], 
        ylim = c(0, 0.08),
        las=2, 
        cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05,col="red")
abline(h= 0.01, col = "blue")
legend("topleft", legend=levels(factor(targets[, class_label])), fill=pal,
       bg="white")

# Also examine the proportion of failed probes (detP >= probe_detP_max).
failed_proportion <- apply(detP, 2, function(x) {
  # All probes detP == 0 means totally failed array.
  if (sum(x) == 0) {
    return(NA)
  } else {
    failed_probes <- sum(x >= 0.01)
    return(failed_probes / length(x))
  }
})
names(failed_proportion) <- targets$Sample_Name

barplot(failed_proportion, col=pal[factor(targets[, class_label])], 
        ylim = c(0, 0.12),
        las=2, 
        cex.names=0.8, ylab="Proportion of failed probes")
abline(h=0.1,col="red")
abline(h=0.05,col="blue")
abline(h= 0.01, col = "green")
legend("topleft", legend=levels(factor(targets[, class_label])), fill=pal,
       bg="white")


# Output the sample sheet, but with an added column with the failed cpg
# proportion.
targets$Prop_failed_probes <- failed_proportion
write.csv(targets, file = file.path(RESULTS_DIR, 
                                    "Sample_sheet_with_prop_failed_cpgs.csv"))


### Remove poor quality samples based on mean detection pval.
keep <- colMeans(detP) < 0.05
print(paste0("Samples to remove based on mean detP threshold of ",
             sample_mean_detP_max,
             ": ",
             paste0(targets$Sample_Name[!keep],";", colnames(detP)[!keep])))
rgSet1 <- rgSet1[,keep]


targets <- targets[keep,]

# remove poor quality samples from detection p-value table
detP <- detP[,keep]


### Remove poor quality samples based on proportion of failed probes.
keep <- apply(detP, 2, function(x) {
  # All probes detP == 0 means totally failed array.
  if (sum(x) == 0) {
    return(FALSE)
  } else {
    failed_probes <- sum(x >= probe_detP_max)
    failed_sample_prop <- failed_probes / length(x)
    return(ifelse(failed_sample_prop < max_failed_probes, TRUE, FALSE))
  }
})
print(paste0("Samples to remove based on proportion of failed probes over ",
             max_failed_probes,
             ": ",
             paste0(targets$Sample_Name[!keep],";", colnames(detP)[!keep])))
rgSet1 <- rgSet1[,keep]
if (!autopilot) {
  print(rgSet)
}
# remove poor quality samples from targets data
targets <- targets[keep,]
if (!autopilot) {
  print(targets[,1:5])
}
# remove poor quality samples from detection p-value table
detP <- detP[,keep]
if (!autopilot) {
  print(dim(detP))
}



################################################################################
## Step 3: Normalization.
if (norm_method == "Quantile"){
  mSetSq <- preprocessQuantile(rgSet1)
} else if (norm_method == "Funnorm") {
  mSetSq <- preprocessFunnorm(rgSet1)
} else if (norm_method == "Noob") {
  mSetSq <- preprocessNoob(rgSet1, dyeMethod = "single")
  # Convert to RatioSet object.
  mSetSq <- ratioConvert(mSetSq)
  # Convert to GenomicRatioSet object.
  mSetSq <- mapToGenome(mSetSq)
} else if (norm_method == "Illumina") {
  mSetSq <- preprocessIllumina(rgSet1, bg.correct = TRUE, normalize = "controls",
                               reference = 1)
  # Convert to RatioSet object.
  mSetSq <- ratioConvert(mSetSq)
  # Convert to GenomicRatioSet object.
  mSetSq <- mapToGenome(mSetSq)
} else {
  stop("[ERROR] normalization method not correctly specified.")
}

# Save raw and normalized b-values as RData files.
# save(rgSet, file = file.path(RESULTS_DIR, "raw.RData"))
# save(mSetSq, file = file.path(RESULTS_DIR, "normalized.RData"))


################################################################################
## Step 4: Filtering.
# Filtering will be performed to normalized beta-values and M-values matrixes.
bVals_original <- getBeta(mSetSq)
mVals_original <- getM(mSetSq)
dim(bVals_original)

detP <- detP[match(rownames(bVals_original),rownames(detP)),]
if (remove_or_NA_probes == "NA") {
  # Mark all probes that have failed (detection pval > threshold) as NAs.
  failed_cpgs <- which(detP > probe_detP_max, arr.ind = TRUE)
  bVals <- bVals_original
  bVals[failed_cpgs] <- NA
  mVals <- mVals_original
  mVals[failed_cpgs] <- NA
  print(paste("NAs before mark failed probes:",
              sum(is.na(bVals_original))))
  print(paste("NAs after marking failed probes:",
              sum(is.na(bVals))))
} else if (remove_or_NA_probes == "remove") {
  # Remove any probes that have failed (detection pval < threshold) in one or 
  # more samples
  keep <- rowSums(detP < probe_detP_max) == dim(bVals_original)[2]
  if (!autopilot) {
    print(paste("Probes to keep:", sum(keep)))
  }
  bVals <- bVals_original[keep, ]
  mVals <- mVals_original[keep, ]
  if (!autopilot) {
    print("Probes before removing failed CpGs.")
    print(dim(bVals_original)[1])
    print("Probes after removing failed CpGs.")
    print(dim(bVals)[1])
  }
} else {
  stop("[ERROR] remove_or_NA_probes var should be set to either NA or remove")
}

# If your data includes males and females, remove probes on the sex chromosomes
if (remove_sex_chr) {
  keep <- !(rownames(bVals) %in% annot$Name[annot$chr %in% c("chrX","chrY")])
  bVals <- bVals[keep, ]
  mVals <- mVals[keep, ]
}
if (!autopilot) {
  print("Probes after removing sex chr CpGs.")
  print(dim(bVals)[1])
}
# Remove probes with SNPs at CpG site and at a single nucleotide extension.
# All are default options.
if (remove_SNP_probes) {
  # Get snp info for all our probes.
  snp_info <- getSnpInfo(mSetSq)
  # Subset for the ones still present.
  snp_info <- snp_info[rownames(bVals), ]
  bVals <- dropLociWithSnps_custom(bVals, snpAnno = snp_info, 
                                   snps = c("CpG", "SBE"), maf = 0)
  mVals <- dropLociWithSnps_custom(mVals, snpAnno = snp_info, 
                                   snps = c("CpG", "SBE"), maf = 0)
}
if (!autopilot) {
  print("Probes after removing probes with SNPs at CpG site or at a single nucleotide extension.")
  print(dim(bVals)[1])
}
# exclude cross reactive probes 
if (remove_cross_reactive) {
  xReactiveProbes <- read.csv(POLYMORPHICEPIC, stringsAsFactors=FALSE)
  keep <- !(rownames(bVals) %in% xReactiveProbes[, 1])
  bVals <- bVals[keep, ]
  mVals <- mVals[keep, ]
}
if (!autopilot) {
  print("Probes after removing cross-reactive probes.")
  print(dim(bVals)[1])
}

# Put correct sample names.
colnames(mVals_original) <- targets$Sample_Name
colnames(bVals_original) <- targets$Sample_Name
colnames(mVals) <- targets$Sample_Name
colnames(bVals) <- targets$Sample_Name



var_bVals <- apply(bVals, 1, var)
cutoff <- mean(var_bVals)
probes_to_keep <- which(ifelse(var_bVals >= cutoff, TRUE, FALSE))
bVals <- bval[probes_to_keep, ]
dim(bVals)
mVals <- mVals[probes_to_keep, ]

## Heatmap
#bVals_corr <- bVals[sample(which(complete.cases(bVals)), 50000),]
dim(bVals_var)
pdf(file = file.path(RESULTS_DIR, "heatmap_LMX_LMH.pdf"))
heatmap.2(bVals_var[1:1000,], col = greenred(75), trace = "none", labRow = FALSE, cexCol = 0.6, 
          ColSideColors = rep(c("#aac5f2","#0748b0") , times=19),
          margins = c(5, 12.5), key.xlab = paste(expression(Beta), "-value", sep = ""), 
          key.title = "", key.ylab = "")
legend(x=0.8, 0.736, legend=c("LMX", "LMH"), col = c("#aac5f2", "#0748b0"), fill= c("#aac5f2", "#0748b0"),bty = "n" )
dev.off()


## Correlation plot (matrix)
corMatrix <- cor(bVals_var)
col <- colorRampPalette(c("#440154", "#482878", "#3e4a89", "#fde725"))
pdf(file = file.path(RESULTS_DIR, "corrMatrix_LMX_LMH.pdf"), height = 15, width = 10)
corrplot::corrplot(corMatrix, cl.lim = c(0.6, 1),
        is.corr = F, method="color",col =col(200), type="lower", tl.col = "black",
        sig.level = 0.01, insig = "blank", tl.cex = 0.5, addCoef.col = "black", number.cex = 0.5)
dev.off()         
idx_LMH <- grep(pattern = "LMH", colnames(bVals))
idx_LMX <- grep(pattern = "LMX", colnames(bVals))

bVal_LMH <- bVals_var[, idx_LMH]
bVal_LMX <- bVals_var[, idx_LMX]
dim(bVal_LMH)

# Correlation dataframe
corr_df <- data.frame()
for (i in 1:ncol(bVal_LMX)){
  test <- cor.test(bVal_LMH[, i], bVal_LMX[, i])
  corrval <- unname(test$estimate)
  pval <- test$p.value
  new <- cbind(as.numeric(corrval), as.numeric(pval))
  corr_df <- rbind(corr_df, new)
}

pdf(file = file.path(RESULTS_DIR, "corr_LMX_LMH.pdf"))
ggplot(corr_df, aes(x=LMX, y=corr_val)) + geom_col(fill="#c4c708", alpha=0.7)+ theme_bw()+
  theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(labels=gsub("LMX", "", corr_df$LMX)) +
  labs(x="", y="Pearson correlation between LMX and LMH") + geom_text(aes(label=round(corr_val, 2)), vjust=-0.45)
dev.off()
colnames(corr_df) <- c("corr_val", "pvalue")
corr_df$LMX <- colnames(bVal_LMX)
corr_df$LMH <- colnames(bVal_LMH)
corr_df <- corr_df[, c(3, 4, 1, 2)]
write.csv(corr_df, file = file.path(RESULTS_DIR, "correlation_LMX_LMH.csv"))

# Mean correlation CpGs
dat <- as.data.frame(cbind(rowMeans(bVal_LMH), rowMeans(bVal_LMX)))
colnames(dat) <- c("LMH", "LMX")

pdf(file = file.path(RESULTS_DIR, "corr_CpGsLMX_LMH.pdf"))
ggplot(dat, aes(x=LMH, y=LMX)) + geom_point(alpha=0.4, color="#0091ff") + theme_classic()+
  annotate(geom = "text", x=0.8, y=0.1, label=paste("Pearson corr.=", round(cor(dat$LMH, dat$LMX), 4))) +
  labs(x="Mean CpG methylation level in patients (LMH)", y="Mean CpG methylation level in patient-derived xenografts (LMX/PDX)")
dev.off()


## Dendogram
dd <- dist(t(bVals_var))
hc <- hclust(dd)
pdf(file = file.path(RESULTS_DIR, "dendogramLMX_LMH.pdf"))
plot(hc)
dev.off()
