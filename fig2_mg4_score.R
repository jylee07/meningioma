# fig2_mg4_score.R: MG4 Meta-Module Score Calculation using ssGSEA

# Load libraries
suppressPackageStartupMessages({
  library(escape)
  library(SingleCellExperiment)
  library(Seurat)
  library(dittoSeq)
  library(GSEABase)
  library(ggplot2)
})

# Load preprocessed Seurat object with normalized data and metadata
# Replace 'seurat_obj.rds' with actual path
seurat_obj <- readRDS("data/seurat_obj.rds")

# Define MG4 meta-module gene sets
MG4 <- list(
  Immune = c("CD14", "TLR8", "LILRA2", "CSF1", "CD1C", "CD74", "HLA-DBQ1", "CD6", "CD4",
             "LILRA1", "HLA-DPA1", "HLA-DPB1", "HLA-DRA", "HLA-DRB5", "HLA-DRB1", "HLA-DQA1",
             "LYVE1", "CCL21", "CD3E", "CD34"),
  NFwt = c("VEGFC", "SEMA5A", "ESM1", "VEGFA", "SERPINE1"),
  Hypermetabolism = c("MDH1", "ACAD9", "MDH2", "DLD", "DECR1", "CYC1", "ACO1", "ACAT1",
                      "CPT2", "ACOX2", "ACOX3", "ACADS"),
  Proliferative = c("CDC20", "MYC", "E2F8", "CCND1", "CHEK1", "CDK1", "TET1", "MKI67",
                    "CDK4", "TOP2A", "BRCA1", "CDKN2A", "FOXM1", "TP53", "MCM2", "MCM6", "MCM4")
)

# Run ssGSEA scoring
mg_scores <- enrichIt(
  obj = seurat_obj,
  gene.sets = MG4,
  method = "ssGSEA",
  groups = 2000,
  cores = 8,
  ssGSEA.norm = TRUE
)

# Calculate axes and groups
mg_scores$nf_immune <- mg_scores$NFwt - mg_scores$Immune
mg_scores$pro_hyper <- mg_scores$Proliferative - mg_scores$Hypermetabolism
mg_scores$D <- apply(mg_scores, 1, function(x) max(x["Immune"], x["NFwt"]) - max(x["Hypermetabolism"], x["Proliferative"]))
mg_scores$x.axis <- ifelse(mg_scores$D > 0, mg_scores$nf_immune, mg_scores$pro_hyper)
mg_scores$group <- NA
mg_scores$group[mg_scores$D > 0 & mg_scores$x.axis > 0] <- "NF2 WT"
mg_scores$group[mg_scores$D > 0 & mg_scores$x.axis < 0] <- "Immune"
mg_scores$group[mg_scores$D < 0 & mg_scores$x.axis < 0] <- "Hypermetabolism"
mg_scores$group[mg_scores$D < 0 & mg_scores$x.axis > 0] <- "Proliferative"

# Plot
p <- ggplot(mg_scores, aes(x = x.axis, y = D, color = group)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  scale_color_manual(values = c(
    "Immune" = "#DE292D",
    "NF2 WT" = "#2475B7",
    "Hypermetabolism" = "#37B64B",
    "Proliferative" = "#F69124"
  )) +
  labs(x = "Relative meta-module score", y = "Relative meta-module score") +
  xlim(c(-1, 1)) + ylim(c(-1, 1)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_vline(xintercept = 0, col = "grey", linetype = "dotted") +
  geom_hline(yintercept = 0, col = "grey", linetype = "dotted")

print(p)

# Save output scores
write.csv(mg_scores, "outputs/fig2_mg4_score.csv")
