#Let's load the full list of neuronal objects arcross studies
library(Seurat)

sc.list <- readRDS('/.../List_of_Individual_SpeciesStudy_Neuronal_Objects_FrogMouseHumanZebrafish_SCTrasformed_Ready_for_Integration.rds')

frog <- sc.list[[1]]
russ <- sc.list[[2]]
kathe <- sc.list[[3]]
gautier <- sc.list[[4]]
yadav <- sc.list[[5]]
## if you want to add zebrafish as well (see Fig4. L-N), it is stored in the list of objects as [[6]]  
# zebrafish <- sc.list[[6]]


frog.vf <- VariableFeatures(frog)
russ.vf <- VariableFeatures(russ)
kathe.vf <- VariableFeatures(kathe)
gautier.vf <- VariableFeatures(gautier)
yadav.vf <- VariableFeatures(yadav)
# zebrafish.vf <- VariableFeatures(zebrafish)

int.features <- Reduce(intersect, list(frog.vf, russ.vf, kathe.vf, gautier.vf, yadav.vf))
# int.features <- Reduce(intersect, list(frog.vf, russ.vf, kathe.vf, gautier.vf, yadav.vf, zebrafish.vf))

# Combine data
sc.list <- c(sc.list[[1]], sc.list[[2]], sc.list[[3]], sc.list[[4]], sc.list[[5]]) 
# sc.list <- c(sc.list[[1]], sc.list[[2]], sc.list[[3]], sc.list[[4]], sc.list[[5]], sc.list[[6]]) 

# Prepare for integration
sc.list.prep <- PrepSCTIntegration(object.list = sc.list, anchor.features = int.features, assay = 'SCT')

# Find integration anchors, for integration via RPCA change reduction to reduction = 'rpca'
sc.anchors <- FindIntegrationAnchors(object.list = sc.list.prep, normalization.method = "SCT", anchor.features = int.features, reduction = 'cca')

# Integrate data, for integration via RPCA, change dims and npcs to 50
sc.integrated <- IntegrateData(anchorset = sc.anchors, normalization.method = "SCT",  dims=1:150)

# Run PCA, UMAP, and find neighbors
sc.integrated <- RunPCA(sc.integrated, npcs = 150)
sc.integrated <- RunUMAP(sc.integrated, dims = 1:150)
sc.integrated <- FindNeighbors(sc.integrated, dims = 1:150)

# Clean up metadata and Save results
saveRDS(sc.integrated, '/.../Neurons_Integrated_Across_FrogMouseHuman_CCA.rds')
