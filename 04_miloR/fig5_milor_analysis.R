# FIG 5: MiloR Differential Abundance and Marker Analysis

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(scater)
  library(scran)
  library(miloR)
  library(tidyverse)
  library(patchwork)
  library(ggrastr)
})

# Convert Seurat to SCE
sce <- as.SingleCellExperiment(lymphocyte)

# Highly variable genes
dec <- modelGeneVar(sce)
hvgs <- getTopHVGs(dec, n = 3000)

# Dimensionality reduction
set.seed(42)
sce <- runPCA(sce, subset_row = hvgs, ncomponents = 11)
sce <- runUMAP(sce, dimred = "PCA", ncomponents = 2)

# Create Milo object
milo <- Milo(sce)
milo <- buildGraph(milo, d = 11, k = 35)
milo <- makeNhoods(milo, k = 35, d = 11, prop = 0.1, refined = TRUE)
plotNhoodSizeHist(milo, bins = 150)

# Design matrix
milo@colData$for_milo <- paste(sce$sample, sce$celltype2, sep = "_")
design_df <- distinct(as_tibble(colData(milo)[, c("sample", "orig.ident")]))
rownames(design_df) <- design_df$sample

# Count and test
milo <- countCells(milo, samples = "sample", meta.data = colData(milo)[, c("sample", "orig.ident", "patient")])
milo <- calcNhoodDistance(milo, d = 30)
res <- testNhoods(milo, design = ~ orig.ident, design.df = design_df)

# Find markers
milo <- buildNhoodGraph(milo)
res$NhoodGroup <- as.numeric(res$logFC < 0)
marker_res <- findNhoodMarkers(milo, res, overlap = 1, compute.new = TRUE)

# DA visualization
plotNhoodGraphDA(milo, res, alpha = 0.1, size_range = c(2, 6), layout = 'TSNE')

# Annotate and classify
anno_res <- annotateNhoods(milo, res, coldata_col = "mad_hclust5")
anno_res$mad_hclust5 <- ifelse(anno_res$mad_hclust5_fraction < 0.7, "Mixed", anno_res$mad_hclust5)
plotDAbeeswarm(anno_res, group.by = "mad_hclust5") + theme_classic()

# Expression matrix prep
milo <- logNormCounts(milo)
keep <- rowSums(logcounts(milo)) != 0
milo <- milo[keep, ]

# Heatmap genes
dec <- modelGeneVar(milo)
hvgs <- getTopHVGs(dec, n = 2000)
marker_df <- findNhoodGroupMarkers(milo, res, subset.row = hvgs, aggregate.samples = TRUE, sample_col = "sample")
rownames(marker_df) <- marker_df$GeneID

# Save top markers of a group
group2_markers <- marker_df[, c("logFC_2", "adj.P.Val_2")]
colnames(group2_markers) <- c("logFC", "adj.P.Val")

# Nhood expression
milo <- calcNhoodExpression(milo, assay = "logcounts", subset.row = hvgs)

# Custom grouping
grouped <- groupNhoods(milo, res, max.lfc.delta = 5, overlap = 5)
grouped$NhoodGroup <- NA
grouped$NhoodGroup <- ifelse((grouped$SpatialFDR < 0.1) & (grouped$logFC < -8), "54", grouped$NhoodGroup)
grouped$NhoodGroup <- ifelse((grouped$SpatialFDR < 0.1) & (grouped$logFC > 8), "70", grouped$NhoodGroup)

plotDAbeeswarm(grouped, group.by = 'NhoodGroup')
