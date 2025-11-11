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

# Load data
frog <- readRDS('frog_SCTranformed_and_Labelled_for_integration_no_downsampling.rds')
mouse <- readRDS('mouse_SCTranformed_and_Labelled_for_integration_no_downsampling.rds')
ariel_neu <- readRDS('mouseRussLevine_SCTranformed_and_Labelled_for_integration_no_downsampling.rds')
yadav <- readRDS('yadav_SCTranformed_and_Labelled_for_integration_no_downsampling.rds')
gautier <- readRDS('gautier_SCTranformed_and_Labelled_for_integration_no_downsampling.rds')
int.features <- readRDS('features_for_frogmousehuman_integration_when_7k_used_for_individualsctranforms.rds')

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
saveRDS(sc.integrated, 'fmh_integrated_no_downsampling_anchors_based_on_vf_overlap_from_sct_with_7k_19_07_24_5datasets.rds')
