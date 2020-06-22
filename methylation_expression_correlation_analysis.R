####################################################################
### Aina Montalban
### Integration between gene expression and DNA methylation on patient-derived xenografts
### Description:
### Intersection between differentially expressed genes (DEGs) and differentially methylated genes (DMPs)
### Find DMP-DEG pairs correlation 
### For the script you should have 4 files:
      ### a) Differentially methylated genes (DMP)
      ### b) Differentially expressed genes (DEG)
      ### c) B-values from DMP
      ### d) expression set from DEG 
##################################

#== Packages
library(ggplot2)
require(gridExtra)

#== Directories
stem_path <- "/home/aina/Downloads" 

RESULTS_DIR <- file.path(
  stem_path, 
  "results")

#DATA_DIR <- file.path(
# stem_path, 
#"data")

setwd(RESULTS_DIR)

#== Read targets
targets <- read.csv(file.path(RESULTS_DIR, "sample_sheet_PDX.csv"), skip = 7)

#== Load beta-values and M-values
load(file.path(RESULTS_DIR, "RES_EXP/meth.RData"))
dim(bVals)
colnames(bVals)

#== Load expression set (w/ replicate averaging)
load(file.path(RESULTS_DIR, "RES_EXP/esetAvg.RData"))
dim(eset_avg)
colnames(eset_avg)


#== Select patients in both datasets (gene expression and DNA methylation)
names_4_bVals <- gsub("A|B", "", colnames(eset_avg))
names_4_bVals <- paste(names_4_bVals, "LMX", sep="")
table(targets$Sample_Group[targets$Sample_Name %in% names_4_bVals])
bVals_exp <- bVals[,names_4_bVals]
dim(bVals_exp) 
targets_exp <- targets[targets$genealogy..LMX..Xenograft.from.liver.mets. %in% names_4_bVals, ]
dim(targets_exp)


#== Load differentially methylated genes (DMP)
DMP <- read.csv("RES_METH/diff-DMPs.csv")
dim(DMP)

# In case we only want to select the DMP located in promoters
promoters <- c("5'UTR", "1stExon", "TSS200", "TSS1500")

DMP <- DMP[grep(paste(promoters,collapse="|"),
               DMP$UCSC_RefGene_Group),]

dim(DMP)
# DMP <- DMP[gsub("\\;.*", "",DMP$UCSC_RefGene_Group) %in% promoters, ]




#== Load differentially expresed genes (DEG)
DEG <- read.csv("RES_EXP/DEGs.csv")
dim(DEG)


#== Intersection
gene_cpg_hyper <- gsub("\\;.*", "",DMP$UCSC_RefGene_Name[DMP$logFC > 0])
gene_cpg_hypo <- gsub("\\;.*", "",DMP$UCSC_RefGene_Name[DMP$logFC < 0])
gene_up <- DEG$geneSymbol[DEG$logFC > 0]
gene_down <- DEG$geneSymbol[DEG$logFC < 0]

# write.csv(gene_cpg_hyper, file = file.path(RESULTS_DIR, "hyper_gene_diff.txt"), quote = FALSE, row.names = FALSE)
# write.csv(gene_cpg_hypo, file = file.path(RESULTS_DIR, "hypo_gene_diff.txt"), quote = FALSE, row.names = FALSE)

length(intersect(gene_cpg_hyper, gene_up))
length(intersect(gene_cpg_hyper, gene_down))
length(intersect(gene_cpg_hypo, gene_down))
length(intersect(gene_cpg_hypo, gene_up))

# Remove DEGs with no genes associated
DEG <- DEG[!is.na(DEG$geneSymbol), ]


# Create an empty data frame, which we will introduce the DEG-DMP pair
DEG_DMP_df <- data.frame()

# Fill data frame with the IDs
for (id in DEG$ID){
  gene_name<- as.character(DEG$geneSymbol[which(DEG$ID == id)]) # select gene name of gene expression dataset
  cpg_names <- as.character(DMP$Name[ grep(gene_name, DMP$UCSC_RefGene_Name)]) # find cpgs related with this name
    if (length(cpg_names) >= 1){ # if there is one or more cpgs associated with the genes
      num <- length(cpg_names) # how may cpgs associated
      new <- cbind(rep(id, num),cpg_names) # add new cpg-gene pairs
      DEG_DMP_df <- rbind(DEG_DMP_df, new)}
}

# Add new colunmns
r_measure <- rep(NA, nrow(DEG_DMP_df)) # correlation value
pvalue <- rep(NA, nrow(DEG_DMP_df)) # p-valie
geneSymbol <- rep(NA, nrow(DEG_DMP_df)) # gene symbol (eg. FOXD3)
rel <- rep(NA, nrow(DEG_DMP_df)) # hyper-down/hyper-up/hypo-up/hypo-down
diff <- rep(NA, nrow(DEG_DMP_df)) # diff. in methylation between RES vs NRES (B > 0.2)
logFC <- rep(NA, nrow(DEG_DMP_df)) # logFC of DEGs


# Loop to calculate the correlation test
for (i in 1: nrow(DEG_DMP_df)){
  pair <- DEG_DMP_df[i, ] # select a pair
  ilm <- as.character(unlist(pair[1])) # save Illumina ID
  cpg <- as.character(unlist(pair[2])) # save cpg id
  test <- cor.test(bVals_exp[cpg,], eset_avg[ilm, ]) # do correlation test between the b-values and log2(mRNA)
  r_measure[i] <- test$estimate # select the corr value
  pvalue[i] <- test$p.value # select the p-value
  
  # Specify more variables (i.e gene symbol, diff. meth, logFC)
  geneSymbol[i] <- as.character(DEG$geneSymbol[which(DEG$ID == ilm)])
  diff[i] <- DMP$Diff_Meth[DMP$Name == cpg]
  logFC[i] <- DEG$logFC[DEG$ID == ilm]
  if (ilm %in% DEG$ID[DEG$logFC > 0] & cpg %in% DMP$Name[DMP$logFC > 0]){
    rel[i] <- "hyper-upregulated"
  } else if (ilm %in% DEG$ID[DEG$logFC < 0] & cpg %in% DMP$Name[DMP$logFC > 0]){
    rel[i] <- "hyper-downregulated"
  } else if (ilm %in% DEG$ID[DEG$logFC > 0] & cpg %in% DMP$Name[DMP$logFC < 0]){
    rel[i] <- "hypo-upregulated"
  } else if(ilm %in% DEG$ID[DEG$logFC < 0] & cpg %in% DMP$Name[DMP$logFC < 0]){
    rel[i] <- "hypo-downregulated"
  }
}

# Add the vectors to the DF
DEG_DMP_df$geneSymbol <- geneSymbol
DEG_DMP_df$corr <- r_measure
DEG_DMP_df$p.value <- pvalue
DEG_DMP_df$rel <- rel
DEG_DMP_df$Diff <- diff
DEG_DMP_df$logFC <- logFC

head(DEG_DMP_df)
dim(DEG_DMP_df)


# Select |r| >= 0.3
DEG_DMP_df_corr <- DEG_DMP_df[abs(DEG_DMP_df$corr) >= 0.3, ]
dim(DEG_DMP_df_corr)

# reorder dataframe
DEG_DMP_df_corr <- DEG_DMP_df_corr[order(DEG_DMP_df_corr$corr, decreasing = TRUE),]
write.csv(DEG_DMP_df_corr, file = "DEG_DMP_corr.csv")

dim(DEG_DMP_df_corr) # check dimensions
table(DEG_DMP_df_corr$rel) # check number of correlations types

head(DEG_DMP_df_corr) # positively correlated
tail(DEG_DMP_df_corr) # negatively correlated




# Plot the 4 positive correlated
plots_pos_corr <- list()
for (i in 1:6){
  pair <- DEG_DMP_df_corr[i, ]
  ilm <- as.character(unlist(pair[1]))
  cpg <- as.character(unlist(pair[2]))
  corval <- as.numeric(unlist(pair[4]))
  geneSymb <- DEG$geneSymbol[which(DEG$ID == ilm)]
  tmp <- as.data.frame(cbind(bVals_exp[cpg,], eset_avg[ilm, ]))
  p <- ggplot(tmp,aes(x=V1, y=V2)) + geom_point() + geom_smooth(method = "lm") +
    labs(x=cpg, y=geneSymb, subtitle = paste(paste("r =", round(corval, 4)))) +
    theme_bw() + theme(plot.subtitle = element_text( hjust = 1)) + theme_classic()
  plots_pos_corr[[i]] <- p
}
png(filename = file.path(RESULTS_DIR, "plots_pos_corr.png"), width = 800, height = 600)
grid.arrange(grobs=plots_pos_corr, ncol=3)
dev.off()



# Plot the negative correlated
plots_neg_corr <- list()
for (j in 0:8){
  pair <- DEG_DMP_df_corr[nrow(DEG_DMP_df_corr) - j, ]
  ilm <- as.character(unlist(pair[1]))
  cpg <- as.character(unlist(pair[2]))
  corval <- as.numeric(unlist(pair[4]))
  geneSymb <- DEG$geneSymbol[which(DEG$ID == ilm)]
  tmp <- as.data.frame(cbind(bVals_exp[cpg,], eset_avg[ilm, ]))
  p <- ggplot(tmp,aes(x=V1, y=V2)) + geom_point() + geom_smooth(method = "lm") +
    labs(x=cpg, y=geneSymb, subtitle = paste(paste("r =", round(corval, 4)))) +
    theme_bw() + theme(plot.subtitle = element_text( hjust = 1)) + theme_classic()
  plots_neg_corr[[j+1]] <- p
}
png(filename = file.path(RESULTS_DIR, "plots_neg_corr_promoters.png"), width = 800, height = 600)
grid.arrange(grobs=plots_neg_corr, ncol=3)
dev.off()


## Gene enrichment analysis
library(enrichR)
results_enrichr <- enrichr(DEG_DMP_df_corr$geneSymbol, databases = c("GO_Biological_Process_2018", "GO_Cellular_Component_2018", "GO_Molecular_Function_2018"))
table_bp <- results_enrichr$GO_Biological_Process_2018
table_cc <- results_enrichr$GO_Cellular_Component_2018
table_mf <- results_enrichr$GO_Molecular_Function_2018
colnames(table_bp)
table_bp <- table_bp[order(table_bp$Adjusted.P.value),]
table_cc <- table_cc[order(table_cc$P.value),]
table_mf <- table_mf[order(table_mf$P.value),]

head(table_bp)
head(table_cc)
head(table_mf)

write.csv(table_bp[table_bp$Adjusted.P.value <= 0.05,], file = "table_bp.csv")
write.csv(table_cc[table_cc$Adjusted.P.value <= 0.05,], file = "table_cc.csv")
write.csv(table_mf[table_mf$Adjusted.P.value <= 0.05,], file = "table_mf.csv")

