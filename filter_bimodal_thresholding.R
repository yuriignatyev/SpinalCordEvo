# Bimodal filtering of low-quality cells based on library complexity
# Manual parameter values (lambda, mu, sigma) set based on visual inspection of log-transformed distributions

library(Seurat)
library(mixtools)

# --- nFeature_RNA filtering ---

x_feat <- log(seu$nFeature_RNA)
x_feat <- x_feat[is.finite(x_feat)]

# Manually chosen initial parameters based on distribution shape
mix_feat <- normalmixEM(
  x = x_feat,
  lambda = c(3, 8),
  mu = c(1, 9),
  sigma = c(1, 1),
  k = 2
)

# Find threshold 
feat_thresh_log <- min(mix_feat$x[which(mix_feat$posterior[, 1] < 0.5)])
feat_thresh <- exp(feat_thresh_log)

plot(mix_feat, which = 2, breaks = 100)
abline(v = feat_thresh_log, col = "blue")

# --- nCount_RNA filtering ---

x_count <- log(seu$nCount_RNA)
x_count <- x_count[is.finite(x_count)]

mix_count <- normalmixEM(
  x = x_count,
  lambda = c(5, 9),
  mu = c(2, 12),
  sigma = c(1, 1),
  k = 2
)

count_thresh_log <- min(mix_count$x[which(mix_count$posterior[, 1] < 0.5)])
count_thresh <- exp(count_thresh_log)

plot(mix_count, which = 2, breaks = 100)
abline(v = count_thresh_log, col = "blue")

# --- subset Seurat object using thresholds ---

seu <- subset(seu, subset = nFeature_RNA > feat_thresh & nCount_RNA > count_thresh)

# --- Seurat normalization and preprocessing ---

seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 7000, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 50)
seu <- RunUMAP(seu, dims = 1:50, verbose = FALSE)
