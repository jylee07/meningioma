## FIG 3: Risk Score Estimation using ssGSEA

# Load libraries
suppressPackageStartupMessages({
  library(escape)
  library(Seurat)
  library(SingleCellExperiment)
  library(GSEABase)
  library(dplyr)
})

# Load data (assumes preprocessed Seurat object or expression matrix named `data`)
# Replace `data` below with the actual Seurat object or SCE expression matrix
# e.g., data <- GetAssayData(seurat_obj, slot = "data")

# Define risk score gene sets
risk_score <- list(
  UP = c("CCL21", "CDC20", "CDKN2A", "CHEK1", "CKS2", "COL1A1", "ESR1", "EZH2", "FBLIM1",
         "FGFR4", "IGF2", "KIF20A", "KRT14", "MDM4", "MMP9", "MYBL1", "PGK1", "PIM1",
         "USF1", "COL6A3", "FOXM1", "COL6A1"),
  DOWN = c("ARID1B", "CCN1", "CCND2", "CD3E", "CDK6", "CDKN2C", "GAS1", "IFNGR1",
           "KDR", "LINC02593", "MUTYH", "PGR", "SPOP", "TAGLN", "TMEM30B")
)

# Run ssGSEA scoring
risk_scores <- enrichIt(
  obj = data,
  gene.sets = risk_score,
  method = "ssGSEA",
  groups = 2000,
  cores = 16,
  ssGSEA.norm = TRUE
)

# Compute difference score (UP - DOWN)
risk_scores$risk_score <- risk_scores$UP - risk_scores$DOWN

# Save result
if (!dir.exists("outputs")) dir.create("outputs")
write.csv(risk_scores, file = "outputs/fig3_risk_score.csv")


