### Delile data (developmental mouse) — counts + metadata → neurons (cardinal types)

#Delile data - All data downloaded from github https://github.com/juliendelile/MouseSpinalCordAtlas/tree/master/dataset
#their expanded metadata from their analysis https://github.com/juliendelile/MouseSpinalCordAtlas/tree/master/output

## paths (edit these)
counts_tsv   <- "/.../UMI_count.tsv"
metadata_tsv <- ".../phenoData_annotated.csv"

## pkgs
library(Seurat)
library(ggplot2)
library(dplyr)
library(tibble)
library(biomaRt)

set.seed(1234)

# ----------load----------
counts <- read.table(counts_tsv, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
colnames(counts) <- gsub("\\.", "-", colnames(counts))  ## match metadata cell IDs

## ---------- Ensembl → gene symbols (keep ambiguous as Ensembl) ----------
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://useast.ensembl.org")
gene_names_ids <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = rownames(counts),
  mart = ensembl
)

duplicates <- gene_names_ids %>%
  group_by(external_gene_name) %>%
  summarise(ensembl_count = n(), .groups = "drop") %>%
  filter(ensembl_count > 1)

exclude_ensembl_ids <- unique(
  gene_names_ids %>%
    filter(external_gene_name %in% duplicates$external_gene_name) %>%
    pull(ensembl_gene_id)
)

mapping <- gene_names_ids %>%
  filter(!(ensembl_gene_id %in% exclude_ensembl_ids)) %>%
  distinct(ensembl_gene_id, external_gene_name) %>%
  deframe()

new_rownames <- ifelse(rownames(counts) %in% names(mapping),
                       mapping[rownames(counts)],
                       rownames(counts))
rownames(counts) <- make.unique(new_rownames)

## ---------- metadata ----------
## note: their file is actually tab-delimited
metadata <- read.delim(metadata_tsv, header = TRUE, quote = "\"", check.names = FALSE, stringsAsFactors = FALSE)
colnames(metadata)[1] <- "cell.id"
rownames(metadata) <- metadata$cell.id

## keep overlap
common_cells <- intersect(colnames(counts), rownames(metadata))
counts   <- counts[, common_cells, drop = FALSE]
metadata <- metadata[common_cells, , drop = FALSE]

## ---------- Seurat object ----------
seu <- CreateSeuratObject(counts = counts, meta.data = metadata)

## ---------- analysis: all cells → neurons ----------
seu <- NormalizeData(seu, verbose = FALSE)
seu <- subset(seu, subset = nFeature_RNA > 1000)
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst")
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 50)

## first 4 PCs = processes → use dims = 5:30
seu <- RunUMAP(seu, dims = 5:30)
seu <- FindNeighbors(seu, dims = 5:30)
seu <- FindClusters(seu, resolution = 0.4)

## quick look (optional)
# DimPlot(seu, label = TRUE, group.by = "replicate_id", label.size = 6, repel = TRUE)
# DimPlot(seu, label = TRUE, label.size = 6, repel = TRUE)

## choose neural clusters (your selection)
neu <- subset(seu, idents = c(1,2,3,7,12,13,14,15,18))
## 1,2,3,14,15,18 = INs; 13 = MNs; 7 = mixed prog/diff/neurons; 12 = diff

## light reprocessing
neu <- NormalizeData(neu)
neu <- FindVariableFeatures(neu, selection.method = "vst", nfeatures = 7000)
neu <- ScaleData(neu, features = rownames(neu))
neu <- RunPCA(neu, npcs = 50, features = VariableFeatures(neu))

## PC1 ribosomal → use dims = 2:30
neu <- RunUMAP(neu, dims = 2:30)
neu <- FindNeighbors(neu, dims = 2:30)
neu <- FindClusters(neu, resolution = 0.9)

## check by replicate (optional)
# DimPlot(neu, label = TRUE, group.by = "replicate_id")

## ---------- integration by replicate (SCT) ----------
sc.list <- SplitObject(neu, split.by = "replicate_id")
sc.list <- lapply(sc.list, function(x) SCTransform(x, variable.features.n = 7000))
features <- SelectIntegrationFeatures(object.list = sc.list, nfeatures = 7000)
sc.list <- PrepSCTIntegration(object.list = sc.list, anchor.features = features)
sc.anchors <- FindIntegrationAnchors(object.list = sc.list, normalization.method = "SCT",
                                     anchor.features = features)
mouse_neu <- IntegrateData(anchorset = sc.anchors, normalization.method = "SCT")

mouse_neu <- RunPCA(mouse_neu, npcs = 50)
mouse_neu <- RunUMAP(mouse_neu, dims = 5:50)
# DimPlot(mouse_neu, label = FALSE, group.by = "replicate_id")

mouse_neu <- FindNeighbors(mouse_neu, dims = 5:50)
mouse_neu <- FindClusters(mouse_neu, resolution = 2)

## remove progenitors (cluster 28) if present
if ("28" %in% levels(mouse_neu)) {
  clusters_keep <- c(levels(mouse_neu)[1:28], levels(mouse_neu)[30:46])
  mouse_neu <- subset(mouse_neu, idents = clusters_keep)
  mouse_neu <- RunPCA(mouse_neu, npcs = 50)
  mouse_neu <- RunUMAP(mouse_neu, dims = 5:50)
  mouse_neu <- FindNeighbors(mouse_neu, dims = 5:50)
  mouse_neu <- FindClusters(mouse_neu, resolution = 2)
}

## rename clusters → neural types (your mapping)
mouse_neu <- RenameIdents(mouse_neu,
  '0'='dI5','1'='dI5','2'='dI5','3'='dI4','4'='dI4','5'='dI4','6'='dI5',
  '7'='7','8'='Differentiating','9'='dI5','10'='V1','11'='dI3','12'='Differentiating','13'='dI1','14'='dI4',
  '15'='Differentiating','16'='V2a','17'='MNs','18'='dI6','19'='Differentiating','20'='Differentiating','21'='dI4',
  '22'='MNs','23'='dI5','24'='V2b','25'='dI5','26'='dI4','27'='V3','28'='dI2','29'='MNs','30'='dI1',
  '31'='Differentiating','32'='V0','33'='dI4','34'='V1','35'='dI5','36'='dI4','37'='dI2','38'='Differentiating',
  '39'='V1','40'='V3','41'='dI4','42'='MNs','43'='dI1','44'='dI5'
)

## fix mixed cluster 7 (V0 vs V1) by subclustering
mouse_neu$neural_types <- Idents(mouse_neu)
Idents(mouse_neu) <- mouse_neu$integrated_snn_res.2  ## adjust if different in your object
if ("7" %in% levels(mouse_neu)) {
  seven <- subset(mouse_neu, idents = "7")
  seven <- RunPCA(seven, npcs = 10)
  seven <- RunUMAP(seven, dims = 2:10)
  seven <- FindNeighbors(seven, dims = 2:10)
  seven <- FindClusters(seven, resolution = 0.2)
  ## rule: subcluster 0 = V0; others = V1
  map7 <- setNames(rep("V1", length(levels(seven))), levels(seven))
  if ("0" %in% names(map7)) map7["0"] <- "V0"
  seven <- RenameIdents(seven, !!!map7)
  mouse_neu$neural_types[Cells(seven)] <- paste(Idents(seven))
  Idents(mouse_neu) <- mouse_neu$neural_types
}

## final levels
levels(mouse_neu) <- c('Differentiating','dI1','dI2','dI3','dI5','dI4','dI6','V0','V1','V2a','V2b','V3','MNs')

## isolate cardinal classes and embed
mouse_neu_ready <- subset(mouse_neu, idents = levels(mouse_neu)[2:13])
mouse_neu_ready <- RunPCA(mouse_neu_ready, npcs = 50)
mouse_neu_ready <- RunUMAP(mouse_neu_ready, dims = 3:50)

## ----------save----------

saveRDS(mouse_neu_ready, "Delile2021_MouseDev_Neurons_CardinalClasses_IntegrationReady.rds")
