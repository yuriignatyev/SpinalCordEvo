###integrating  developmental data across species

library(Seurat)
set.seed(54)



## mark species
mouse$species <- 'Mouse'
frog_m$species <- 'Frog'   

## leave only cardinal classes
cardinal_classes <- c('dI1','dI2','dI3','dI4','dI5','dI6','V0','V1','V2a','V2b','V3','MNs')

mouse <- readRDS(/.../Delile2021_MouseDev_Neurons_CardinalClasses_IntegrationReady.rds)
frog_m <- readRDS(/.../Xenopus54_Neurons_CardinalClasses_MouseOrthologs_IntegrationReady.rds)

mouse  <- subset(mouse,  idents = cardinal_classes)
frog_m <- subset(frog_m, idents = cardinal_classes)

## integration
sc.list <- list(frog_m, mouse)   

sc.list <- lapply(sc.list, function(x) {
  SCTransform(x, variable.features.n = 7000, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = sc.list, nfeatures = 7000)

## remove Rgs/Rpl genes from features 
features <- setdiff(features, grep("^(Rgs|Rpl)", features, value = TRUE)) 

sc.list <- PrepSCTIntegration(object.list = sc.list, anchor.features = features, verbose = FALSE)

sc.anchors <- FindIntegrationAnchors(object.list = sc.list,
                                     normalization.method = "SCT",
                                     anchor.features = features,
                                     verbose = FALSE)

scint <- IntegrateData(anchorset = sc.anchors, normalization.method = "SCT", verbose = FALSE)
scint <- RunPCA(scint, npcs = 50, verbose = FALSE)
scint <- RunUMAP(scint, dims = 2:50, verbose = FALSE)

saveRDS(scint, "Xenopus_Mouse_Dev_Integration_CCA_Main_Cardinal_Classes.rds")

scint <- FindNeighbors(scint, dims = 2:50, verbose = FALSE)
scint <- FindClusters(scint, resolution = 0.9, verbose = FALSE)


