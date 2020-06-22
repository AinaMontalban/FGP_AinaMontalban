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

library(ggplot2)
library(enrichR)
library(stringr)
#Define directories
stem_path <- "/home/amontalban/Downloads" 

# B-values stored in the results directory
RESULTS_DIR <- file.path(stem_path, "results")

# Path where project data (as samplesheet.csv) is stored.
DATA_DIR <- file.path(stem_path, "results")
EXP_DIR <- file.path(stem_path, "RES_EXP")
METH_DIR <- file.path(stem_path, "RES_METH")

# GO terms from DNA methylation signatures
table_bp_meth <- read.csv(file.path(METH_DIR, "GO_Biological_Process.csv"))
table_cc_meth <- read.csv(file.path(METH_DIR, "GO_Cellular_Component.csv"))
table_mf_meth <- read.csv(file.path(METH_DIR, "GO_Molecular_Function.csv"))

table_bp_meth <- table_bp_meth[table_bp_meth$Adjusted.P.value <= 0.05,]
table_mf_meth <- table_mf_meth[table_mf_meth$Adjusted.P.value <= 0.05,]
table_cc_meth <- table_cc_meth[table_cc_meth$Adjusted.P.value <= 0.05,]
table_bp_meth$Group <- "BP"
table_cc_meth$Group <- "CC"
table_mf_meth$Group <- "MF"
GO_meth <- rbind(table_bp_meth, table_cc_meth, table_mf_meth)

ggplot(GO_meth, aes(y=as.numeric(gsub("\\/.*", "", Overlap)), x=reorder(Term, -Adjusted.P.value), fill=Group)) + 
  geom_col(width = 0.9, show.legend = F) +  labs(y="Number of genes", x="", caption = "") + theme_classic() +
  geom_text(aes(label=gsub("\\/.*", "", Overlap)), position="stack", hjust=-0.3, color="grey", size=3) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,350)) + scale_fill_brewer(palette = 3) +
  theme(axis.line.x = element_blank(), axis.text.x = element_text(size = 3.6),panel.grid = element_blank(),
        strip.background = element_rect(color = "white")) + 
  facet_grid(Group~., scale="free_y", space = "free") + coord_flip()


# Hypermethylated
table_bp_hyper <- read.csv(file.path(METH_DIR, "GO_Biological_Process_hyper.csv"))
table_cc_hyper <- read.csv(file.path(METH_DIR, "GO_Cellular_Component_hyper.csv"))
table_mf_hyper <- read.csv(file.path(METH_DIR, "GO_Molecular_Function_hyper.csv"))

table_bp_hyper <- table_bp_hyper[table_bp_hyper$Adjusted.P.value <= 0.05,]
table_mf_hyper <- table_mf_hyper[table_mf_hyper$Adjusted.P.value <= 0.05,]
table_cc_hyper <- table_cc_hyper[table_cc_hyper$Adjusted.P.value <= 0.05,]

# Hypomethylated
table_bp_hypo <- read.csv(file.path(METH_DIR, "GO_Biological_Process_hypo.csv"))
table_cc_hypo <- read.csv(file.path(METH_DIR, "GO_Cellular_Component_hypo.csv"))
table_mf_hypo <- read.csv(file.path(METH_DIR, "GO_Molecular_Function_hypo.csv"))

table_bp_hypo <- table_bp_hypo[table_bp_hypo$Adjusted.P.value <= 0.05,]
table_mf_hypo <- table_mf_hypo[table_mf_hypo$Adjusted.P.value <= 0.05,]
table_cc_hypo <- table_cc_hypo[table_cc_hypo$Adjusted.P.value <= 0.05,]

table_bp_hyper$Group <- "HYPER"
table_bp_hypo$Group <- "HYPO"
GOS <- rbind(table_bp_hypo[1:5,], table_bp_hyper[1:5,])

ggplot(GOS, aes(y=as.numeric(Combined.Score), x=reorder(Term, -Adjusted.P.value), fill=Group)) + 
  geom_col(width = 0.9) +  labs(y="Combined Score", x="", subtitle = "") + theme_classic() +
  geom_text(aes(label=gsub("\\/.*", "", Overlap)), position="stack", hjust=-0.3, color="grey", size=3) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,60)) + scale_fill_manual(values = c("#cc0000", "#2d8658"), name="",
                                                                            labels=c("Hypermethylated genes in NRES", "Hypomethylated genes in NRES")) +
  theme(axis.line.x = element_blank(), axis.text.x = element_text(size = 6),panel.grid = element_blank(),
        strip.background = element_rect(color = "white"), axis.text.y = element_text(size = 11.5), 
        plot.subtitle = element_text(hjust = 1.07, face = "bold"), strip.text = element_blank(), legend.position = "top") + 
  facet_grid(Group~., scale="free_y", space = "free", 
             labeller = labeller(Group=list("Biological Process", "Cellular Component", "Molecular FUnction"))) +
  coord_flip()


# transcriptome signatues
# up-regulated
table_bp_up <- read.csv(file.path(EXP_DIR, "GO_Biological_Process_up.csv"))
table_cc_up <- read.csv(file.path(EXP_DIR, "GO_Cellular_Component_up.csv"))
table_mf_up <- read.csv(file.path(EXP_DIR, "GO_Molecular_Function_up.csv"))

table_bp_up <- table_bp_up[table_bp_up$Adjusted.P.value <= 0.05,]
table_mf_up <- table_mf_up[table_mf_up$Adjusted.P.value <= 0.05,]
table_cc_up <- table_cc_up[table_cc_up$Adjusted.P.value <= 0.05,]


# down-regulated
table_bp_down <- read.csv(file.path(EXP_DIR, "GO_Biological_Process_down.csv"))
table_cc_down <- read.csv(file.path(EXP_DIR, "GO_Cellular_Component_down.csv"))
table_mf_down <- read.csv(file.path(EXP_DIR, "GO_Molecular_Function_down.csv"))

table_bp_down <- table_bp_down[table_bp_down$Adjusted.P.value <= 0.05,]
table_mf_down <- table_mf_down[table_mf_down$Adjusted.P.value <= 0.05,]
table_cc_down <- table_cc_down[table_cc_down$Adjusted.P.value <= 0.05,]
table_bp_up$Group <- "UP"
table_bp_down$Group <- "DOWN"

GOS_exp <- rbind(table_bp_up[1:5, ],table_bp_down)
ggplot(GOS_exp, aes(y=as.numeric(Combined.Score), x=reorder(Term, -Adjusted.P.value), fill=Group)) + 
  geom_col(width = 0.9) +  labs(y="Combined Score", x="", subtitle = "") + theme_classic() +
  geom_text(aes(label=gsub("\\/.*", "", Overlap)), position="stack", hjust=-0.3, color="grey", size=3) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) + scale_fill_brewer(palette = 1, name="",
                                                                             labels=c("Down-regulated genes in NRES", "Up-regulated genes in NRES")) +
  theme(axis.line.x = element_blank(), axis.text.x = element_text(size = 6),panel.grid = element_blank(),
        strip.background = element_rect(color = "white"), axis.text.y = element_text(size = 11.5), 
        plot.subtitle = element_text(hjust = 1.07, face = "bold"), strip.text = element_blank(), legend.position = "top") + 
  facet_grid(Group~., scale="free_y", space = "free", 
             labeller = labeller(Group=list("Biological Process", "Cellular Component", "Molecular FUnction"))) +
  coord_flip()
