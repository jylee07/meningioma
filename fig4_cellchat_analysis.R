# fig4_cellchat_primary_vs_recurrent.R

# Load required libraries
suppressPackageStartupMessages({
  library(CellChat)
  library(patchwork)
  library(Seurat)
  library(SeuratObject)
})

options(stringsAsFactors = FALSE)
future::plan("multisession", workers = 4)

# Define helper function
run_cellchat_pipeline <- function(seurat_obj) {
  expr_data <- GetAssayData(seurat_obj, slot = "data")
  meta_data <- data.frame(labels = seurat_obj$celltype2, row.names = names(Idents(seurat_obj)))
  cellchat <- createCellChat(object = expr_data, meta = meta_data, group.by = "labels", datatype = "RNA")
  cellchat <- addMeta(cellchat, meta = meta_data)
  cellchat <- setIdent(cellchat, ident.use = "labels")
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  return(cellchat)
}

# Load preprocessed Seurat objects (tumor_primary and tumor_recurrent)
# These should be loaded prior to this script

# Run CellChat
cellchat_primary <- run_cellchat_pipeline(tumor_primary)
cellchat_recurrent <- run_cellchat_pipeline(tumor_recurrent)

# Update cluster label order
new_order <- c(
  "Active monocyte", "C1Q+ macrophage", "Dendritic cell", "M1 macrophage", "M2 macrophage",
  "Macrophage", "Monocyte", "Immune", "NF2 wt", "Hypermetabolism", "Proliferative", "Pro macrophage"
)
cellchat_primary <- updateClusterLabels(cellchat_primary, new.order = new_order)
cellchat_recurrent <- updateClusterLabels(cellchat_recurrent, new.order = new_order)

# Visualization - Bubble plot
pdf("results/fig4_primary_tumor_bubble.pdf", width = 20, height = 16)
netVisual_bubble(cellchat_primary, sources.use = 8:11, targets.use = c(1:7, 12),
                 remove.isolate = FALSE, sort.by.source = TRUE) +
  theme(axis.text.x = element_text(angle = 45, face = "bold", size = 11, vjust = 1))
dev.off()

pdf("results/fig4_recurrent_tumor_bubble.pdf", width = 20, height = 16)
netVisual_bubble(cellchat_recurrent, sources.use = 8:11, targets.use = c(1:7, 12),
                 remove.isolate = FALSE, sort.by.source = TRUE) +
  theme(axis.text.x = element_text(angle = 45, face = "bold", size = 11, vjust = 1))
dev.off()

# Merge and compare
cellchat_merged <- mergeCellChat(list(Primary = cellchat_primary, Recurrent = cellchat_recurrent), add.names = TRUE)

# Functional similarity
cellchat_merged <- computeNetSimilarityPairwise(cellchat_merged, type = "functional")
cellchat_merged <- netEmbedding(cellchat_merged, type = "functional")
cellchat_merged <- netClustering(cellchat_merged, type = "functional")
netVisual_embeddingPairwise(cellchat_merged, type = "functional", label.size = 3.5)

# Structural similarity
cellchat_merged <- computeNetSimilarityPairwise(cellchat_merged, type = "structural")
cellchat_merged <- netEmbedding(cellchat_merged, type = "structural")
cellchat_merged <- netClustering(cellchat_merged, type = "structural")
netVisual_embeddingPairwise(cellchat_merged, type = "structural", label.size = 3.5)
netVisual_embeddingPairwiseZoomIn(cellchat_merged, type = "structural", nCol = 2)

# Rank and compare
rankNet(cellchat_recurrent, mode = "single", measure = "weight",
        sources.use = 4, targets.use = c(1:3, 6:9, 11), stacked = FALSE, do.stat = TRUE)

rankSimilarity(cellchat_merged, type = "functional")

gg1 <- rankNet(cellchat_merged, mode = "comparison", measure = "weight",
               sources.use = 12, targets.use = c(1:3, 6:9, 11), stacked = TRUE, do.stat = TRUE)
gg2 <- rankNet(cellchat_merged, mode = "comparison", measure = "weight",
               sources.use = 12, targets.use = c(1:3, 6:9, 11), stacked = FALSE, do.stat = TRUE)

gg1 + gg2
