# FIG 2: Pseudo-bulk DEG and pathway enrichment analysis

# Load libraries
library(Seurat)
library(SingleCellExperiment)
library(scater)
library(Matrix.utils)
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(apeglm)

# Prepare input
counts <- cancer@assays$RNA@counts
metadata <- cancer@meta.data
sce <- SingleCellExperiment(assays = list(counts = counts), colData = metadata)

# QC
sce_qc <- perCellQCMetrics(sce)
sce$is_outlier <- isOutlier(sce_qc$total, nmads = 2, type = "both", log = TRUE)
sce <- sce[, !sce$is_outlier]
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]

# Aggregate counts
groups <- colData(sce)[, c("sample", "orig.ident")]
pb <- aggregate.Matrix(t(counts(sce)), groupings = groups, fun = "sum")
pb <- t(pb)

# Metadata
sample_ids <- colnames(pb)
group_labels <- rep(c("primary", "recurrent"), each = 7)
meta <- data.frame(sample = sample_ids, orig.ident = group_labels)
rownames(meta) <- meta$sample

# DESeq2
dds <- DESeqDataSetFromMatrix(countData = pb, colData = meta, design = ~ orig.ident)
dds <- DESeq(dds)
rld <- rlog(dds, blind = TRUE)
plotPCA(rld, intgroup = "orig.ident")
pheatmap(cor(assay(rld)), annotation_col = meta[, "orig.ident", drop = FALSE])

# DEG results
res <- results(dds, contrast = c("orig.ident", "recurrent", "primary"), alpha = 0.05)
res <- lfcShrink(dds, contrast = c("orig.ident", "recurrent", "primary"), res = res, type = "normal")
res_tbl <- as_tibble(data.frame(gene = rownames(res), res))

# Enrichment
gmt <- read.gmt("data/c2.all.v7.1.symbols.gmt.txt")
genes_up <- res_tbl %>% filter(padj < 0.05, log2FoldChange > 1) %>% pull(gene)
genes_down <- res_tbl %>% filter(padj < 0.05, log2FoldChange < -1) %>% pull(gene)
recurr_enrich <- enricher(gene = genes_up, TERM2GENE = gmt, pAdjustMethod = "fdr", universe = res_tbl$gene)
primary_enrich <- enricher(gene = genes_down, TERM2GENE = gmt, pAdjustMethod = "fdr", universe = res_tbl$gene)

# Save results
write.csv(res_tbl, "results/fig2_pseudobulk_deg_results.csv", row.names = FALSE)
write.csv(recurr_enrich@result, "results/fig2_recurr_enrichment.csv", row.names = FALSE)
write.csv(primary_enrich@result, "results/fig2_primary_enrichment.csv", row.names = FALSE)
