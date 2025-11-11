#!/usr/bin/env Rscript

# Load required libraries
library(optparse)
library(Seurat)
library(glmGamPoi)

# Define command line arguments
option_list <- list(
  make_option(c("--data_dir"), type="character", default=NULL, help="Path to the data directory", metavar="character")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if mandatory arguments are provided
if (is.null(opt$data_dir)) {
  stop("You must provide a data directory path with --data_dir")
}

# Set working directory
setwd(opt$data_dir)



frog <- readRDS('/nfs/scistore23/sweengrp/iignatev/Tetrapod_SC/Fig4/Adult_Frog_Dataset_Ready_For_Integration_Not_SCTransformed.rds')
mouse <- readRDS('/nfs/scistore23/sweengrp/iignatev/Tetrapod_SC/Fig4/Adult_Mouse_KatheCourtine_Dataset_Ready_For_Integration_Not_SCTransformed.rds')
ariel_neu <- readRDS('/nfs/scistore23/sweengrp/iignatev/Tetrapod_SC/Fig4/Adult_Mouse_RussLevine_Dataset_Ready_For_Integration_Not_SCTransformed.rds')
yadav <- readRDS('/nfs/scistore23/sweengrp/iignatev/Tetrapod_SC/Fig4/Adult_Human_Dataset_Yadav_Ready_For_Integration_Not_SCTransformed.rds')
gautier <- readRDS('/nfs/scistore23/sweengrp/iignatev/Tetrapod_SC/Fig4/Adult_Human_Dataset_Gautier_Ready_For_Integration_Not_SCTransformed.rds')

ariel_neu$species <- 'Mouse'
ariel_neu$coarse <- readRDS('/nfs/scistore23/sweengrp/iignatev/Tetrapod_SC/Fig4/mouse_russ_coarse_ids.rds')

frog <- SCTransform(frog, vars.to.regress = c('nCount_RNA', 'nFeature_RNA'), variable.features.n = 10000,vst.flavor = "v2")
mouse <- SCTransform(mouse, vars.to.regress =  c('nCount_RNA', 'nFeature_RNA', 'batch'), variable.features.n = 10000,vst.flavor = "v2")
ariel_neu <- SCTransform(ariel_neu, vars.to.regress =  c('nCount_RNA', 'nFeature_RNA'), variable.features.n = 10000,vst.flavor = "v2")
yadav <- SCTransform(yadav, vars.to.regress =  c('nCount_RNA', 'nFeature_RNA', 'sample'), variable.features.n = 10000,vst.flavor = "v2")
gautier <- SCTransform(gautier, vars.to.regress =  c('nCount_RNA', 'nFeature_RNA'), variable.features.n = 10000,vst.flavor = "v2")


frog.vf <- VariableFeatures(frog)
mouse.vf <- VariableFeatures(mouse)
ariel_neu.vf <- VariableFeatures(ariel_neu)
yadav.vf <- VariableFeatures(yadav)
gautier.vf <- VariableFeatures(gautier)

int.features <- intersect(frog.vf, mouse.vf)
int.features <- intersect(int.features, ariel_neu.vf)
int.features <- intersect(int.features, yadav.vf)
int.features <- intersect(int.features, gautier.vf)

# Combine data
sc.list_1 <- c(frog, c(mouse, ariel_neu))
sc.list_2 <- c(yadav, gautier)
sc.list <- c(sc.list_1, sc.list_2)

# Prepare for integration
sc.list.prep <- PrepSCTIntegration(object.list = sc.list, anchor.features = int.features, assay = 'SCT')

# Find integration anchors
sc.anchors <- FindIntegrationAnchors(object.list = sc.list.prep, normalization.method = "SCT", anchor.features = int.features, reduction = 'cca')

# Integrate data
sc.integrated <- IntegrateData(anchorset = sc.anchors, normalization.method = "SCT",  dims=1:150)

# Run PCA, UMAP, and find neighbors
sc.integrated <- RunPCA(sc.integrated, npcs = 150)
sc.integrated <- RunUMAP(sc.integrated, dims = 1:150)
sc.integrated <- FindNeighbors(sc.integrated, dims = 1:150)

# Save results
saveRDS(sc.integrated, 'fmh_integrated_no_downsampling_anchors_based_on_vf_overlap_from_sct_with_10k_19_07_24_5datasets_1.rds')
