## Frog adult spinal cord nuclei + spatial ##

library(Seurat)
library(dplyr)
library(ggplot2)


## 1. Load and merge adult frog nuclei    ##

## 3 frog adult spinal cords
frog_081 <- Read10X('/mnt/fas01fs/fas01fs_group/biogrp/Sweeney/DA0016/20210209_15_Xla_Adult_spinalcordnuclei_combined/20210209_15_Xla_Adult_spinalcordnuclei_combined/cellranger-introns/count_146081/outs/filtered_feature_bc_matrix')
frog_126 <- Read10X('/mnt/fas01fs/fas01fs_group/biogrp/Sweeney/DA0016/20220831_xla_adult_spinalcord/cellranger-introns/count_205126/outs/filtered_feature_bc_matrix')
frog_104 <- Read10X('/mnt/fas01fs/fas01fs_group/biogrp/Sweeney/DA0016/20220831_xla_adult_spinalcord/cellranger-introns/count_205104/outs/filtered_feature_bc_matrix')


frog_081 <- ReadMtx(
  mtx        = file.path(data_dir, "146081_matrix.mtx.gz"),
  features   = file.path(data_dir, "146081_features.tsv.gz"),
  cells      = file.path(data_dir, "146081_barcodes.tsv.gz")
)
frog_081 <- CreateSeuratObject(frog_081)

frog_126 <- ReadMtx(
  mtx        = file.path(data_dir, "205126_matrix.mtx.gz"),
  features   = file.path(data_dir, "205126_features.tsv.gz"),
  cells      = file.path(data_dir, "205126_barcodes.tsv.gz")
)
frog_126 <- CreateSeuratObject(frog_126)

frog_104 <- ReadMtx(
  mtx        = file.path(data_dir, "205104_matrix.mtx.gz"),
  features   = file.path(data_dir, "205104_features.tsv.gz"),
  cells      = file.path(data_dir, "205104_barcodes.tsv.gz")
)
frog_104 <- CreateSeuratObject(frog_104)


frog_081$replicate <- '1'
frog_126$replicate <- '2'
frog_104$replicate <- '3'

# combine them 
frog <- merge(frog_081, c(frog_126, frog_104))
frog[["RNA"]] <- JoinLayers(frog[["RNA"]])

## show batch effects for supp fig
frog <- NormalizeData(frog)
frog <- FindVariableFeatures(frog, selection.method = 'vst', nfeatures = 7000)
frog <- ScaleData(frog, features = VariableFeatures(frog))
frog <- RunPCA(frog, npcs = 50)
frog <- RunUMAP(frog, dims = 1:50)

## hist of nFeature_RNA
pdf('nFeature_RNA_by_replicate_densities.pdf', width = 7, height = 5)
ggplot(frog@meta.data, aes(x = nFeature_RNA, fill = replicate)) + 
  geom_density(alpha = 0.3, adjust = 1.5) +
  theme_minimal() +
  theme(
    text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(hjust = 1, size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )
dev.off()


## 2. Integration by replicate (all cells) ##

sc.list <- c(frog_081, c(frog_104, frog_126))

sc.list <- lapply(X = sc.list, FUN = function(x) {
  SCTransform(x, variable.features.n = 7000)
})

features <- SelectIntegrationFeatures(object.list = sc.list, nfeatures = 7000)

## remove ribosomal / mitochondrial / generic metabolic genes from features
features <- features[!grepl("rgs|rpl|COX|ATP|CYT|ND|16S|12S", features)]

sc.list <- PrepSCTIntegration(object.list = sc.list, anchor.features = features)

sc.anchors <- FindIntegrationAnchors(
  object.list = sc.list,
  normalization.method = "SCT",
  anchor.features = features
)


scint <- IntegrateData(anchorset = sc.anchors, normalization.method = "SCT")
scint <- RunPCA(scint, npcs = 50)
scint <- RunUMAP(scint, dims = 1:50)
scint <- FindNeighbors(scint, dims = 1:50)
scint <- FindClusters(scint, resolution = 0.5)

## cluster markers (for annotation)
scint.markers <- FindAllMarkers(scint, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

split_data <- scint.markers %>%
  group_split(cluster)

sorted_data <- lapply(split_data, function(df) {
  df %>% arrange(p_val_adj)
})

## 3. Annotate broad cell types           



scint <- RenameIdents(
  scint,
  '22' = 'MNs', '27' = 'CSF-cNs',
  '0' = 'INs', '7' = 'INs', '8' = 'INs', 
  '9' = 'INs', '12' = 'INs', '16' = 'INs', '18' = 'INs',
  '19' = 'INs', '20' = 'INs', '21' = 'INs', '23' = 'INs',
  '24' = 'INs', '25' = 'INs', '28' = 'INs', '29' = 'INs',
  '30' = 'INs',
  '6' = 'OPCs',
  '4' = 'Astro',
  '1' = 'OCs', '2' = 'OCs', '3' = 'OCs', '13' = 'OCs',
  '11' = 'Endo'
)


##clean up metadata and save
saveRDS(scint, '/.../AllCells_Integrated_Xenopus_Adult.rds')

## 5. Plots for AdultAllCells object       

##rename remaning clusters into 'Other'

## the object is pre-labeled and ordered for plotting

DotPlot(
  scint,
  features = rev(c(
    'elavl4.L', 'nrxn3.L', 'gad2.L', 'slc17a7.S',
    'slit2.S', 'pkd2l1.L', 'plp1.L', 'mbp.L', 'slc4a4.L',
    'etnppl.L', 'pdgfra.L', 'prx.L', 'mpz.L', 'slc2a2.L',
    'klf7.L', 'lyn.L', 'ptprc.L'
  )),
  cols = c('white', 'black'),
  dot.scale = 10
) +
  RotatedAxis() +
  coord_flip()



DimPlot(
  scint,
  pt.size = 0.5,
  reduction = "umap",
  label = TRUE,
  label.size = 8.5,
  cols = c(
    'INs' = 'skyblue2', 'MNs' = 'dodgerblue4',
    'CSF-cNs' = 'cornflowerblue', 'OCs' = 'orange',
    'Astro' = 'salmon', 'OPCs' = 'coral1',
    'Schwann' = 'sienna1', 'Endo' = 'bisque3',
    'Microglia' = 'lavenderblush4', 'Other' = 'gray95'
  ),
  repel = TRUE
) +
  theme(legend.text = element_text(color = "grey3", size = 15)) +
  NoLegend() +
  NoAxes()



##Neurons
## 6. Reintegrate neurons only             


neu <- subset(scint, idents = c(18, 24, 19, 0, 28, 29, 12, 7, 8, 9, 23, 16, 25, 30, 21, 20, 22))
DefaultAssay(neu) <- 'RNA'
Idents(neu) <- neu$replicate

neu_1 <- subset(neu, idents = '1')
neu_2 <- subset(neu, idents = '2')
neu_3 <- subset(neu, idents = '3')

neu_1 <- GetAssayData(neu_1, layer = 'counts.1')
neu_2 <- GetAssayData(neu_2, layer = 'counts.3')
neu_3 <- GetAssayData(neu_3, layer = 'counts.2')

neu_1 <- CreateSeuratObject(neu_1)
neu_2 <- CreateSeuratObject(neu_2)
neu_3 <- CreateSeuratObject(neu_3)

sc.list <- c(neu_1, c(neu_2, neu_3))

sc.list <- lapply(X = sc.list, FUN = function(x) {
  SCTransform(x, variable.features.n = 7000)
})

features <- SelectIntegrationFeatures(object.list = sc.list, nfeatures = 7000)

## remove ribosomal / mitochondrial etc.
features <- features[!grepl("rgs|rpl|COX|ATP|CYT|ND|16S|12S", features)]

sc.list <- PrepSCTIntegration(object.list = sc.list, anchor.features = features)

sc.anchors <- FindIntegrationAnchors(
  object.list = sc.list,
  normalization.method = "SCT",
  anchor.features = features
)

scint <- IntegrateData(anchorset = sc.anchors, normalization.method = "SCT")
scint <- RunPCA(scint, npcs = 50)
scint <- RunUMAP(scint, dims = 1:50)
scint <- FindNeighbors(scint, dims = 1:50)
scint <- FindClusters(scint, resolution = 1.5)

### differential markers
DefaultAssay(scint) <- 'RNA'
scint[["RNA"]] <- JoinLayers(scint[["RNA"]])
scint <- NormalizeData(scint)

neu.markers <- FindAllMarkers(scint, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

neu.markers.split <- neu.markers %>%
  group_split(cluster)

neu.markers.sorted <- lapply(neu.markers.split, function(df) {
  df %>% arrange(p_val_adj)
})


## 7. Remove low-quality neurons  (such as cluster 0)         ##

zero <- subset(scint, idents = 0)
DefaultAssay(zero) <- 'integrated'
zero <- RunPCA(zero, npcs = 10)
zero <- RunUMAP(zero, dims = 1:10)
zero <- FindNeighbors(zero, dims = 1:10)
zero <- FindClusters(zero, resolution = 0.5)

DefaultAssay(zero) <- 'RNA'
zero.markers <- FindAllMarkers(zero, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

zero.markers.split <- zero.markers %>%
  group_split(cluster)

zero.markers <- lapply(zero.markers.split, function(df) {
  df %>% arrange(p_val_adj)
})

zero <- RenameIdents(
  zero,
  '0' = 'low', '1' = 'good', '2' = 'good', '3' = 'good',
  '4' = 'good', '5' = 'good', '6' = 'good', '7' = 'good'
)

scint$quality <- 'good'
scint$quality[Cells(zero)] <- paste(Idents(zero))

DimPlot(scint, label = TRUE, group.by = 'quality')

DefaultAssay(scint) <- 'integrated'
Idents(scint) <- scint$quality

scint <- subset(scint, idents = 'good')
scint <- RunPCA(scint, npcs = 50)
scint <- RunUMAP(scint, dims = 1:50)
scint <- FindNeighbors(scint, dims = 1:50)

##This is where we got coarse neuronal groups
scint <- FindClusters(scint, resolution = 0.07)

#And then start labelling neural types

## 8. Excitatory dorsal populations        ##


exc <- subset(scint, idents = c(8, 6, 3, 5))
exc <- RunPCA(exc, npcs = 50)
exc <- RunUMAP(exc, dims = 1:50)
exc <- FindNeighbors(exc, dims = 1:50)
exc <- FindClusters(exc, resolution = 1)

genes_exc_dors <- rev(c(
  'LOC108713448', 'slc17a7.L', 'lmx1b.1.S', 'cpne4.L',
  'bnc2.S', 'nmu.L.1', 'skor1.L', 'npy.L', 'tac1.S',
  'satb1.L', 'cdh23.L', 'tcf4.S', 'penk.L'
))

DotPlot(exc, features = genes_exc_dors, assay = 'RNA') +
  RotatedAxis() +
  coord_flip() +
  scale_colour_gradient2(low = "white", mid = "gray85", high = "black")

exc <- RenameIdents(
  exc,
  '4' = 'E1', '11' = 'E2', '14' = 'E3', '13' = 'E4', '3' = 'E5', '6' = 'E6',
  '5' = 'E7', '0' = 'E8', '9' = 'E9', '1' = 'E10', '8' = 'E11', '12' = 'E12',
  '2' = 'E13', '7' = 'E14', '10' = 'E15', '15' = 'E16', '16' = 'E17'
)

exc <- RenameIdents(
  exc,
  'E10' = 'E8', 'E11' = 'E9', 'E12' = 'E10',
  'E9' = 'E11', 'E8' = 'E12'
)

levels(exc) <- c(
  'E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E8', 'E9',
  'E10', 'E11', 'E12', 'E13', 'E14', 'E15', 'E16', 'E17'
)


## 9. Ventral / midline excit/inhib        ##


vent <- subset(scint, idents = 0)
vent <- RunPCA(vent, npcs = 50)
vent <- RunUMAP(vent, dims = 1:50)
vent <- FindNeighbors(vent, dims = 1:50)
vent <- FindClusters(vent, resolution = 0.6)

levels(vent) <- c(
  13, 16, 10, 19, 14, 0, 6, 5, 9, 2, 7, 3, 11, 17, 12, 8, 18, 1, 15, 4, 20
)

genes_vent <- c(
  'zfhx3.L', 'LOC108713448', 'slc17a7.L', 'gad2.L', 'slc6a5.S',
  'zfhx4.L', 'lmx1b.1.S', 'pou2f2.L', 'prdm16.S', 'lhx2.S', 'lhx9.S',
  'vsx2.S', 'sim1.S', 'hoxc10.L', 'hoxc9.L', 'esrrg.S', 'foxp2.L',
  'foxp1.S', 'irx5.L', 'nr4a2.L', 'en1.L', 'tfap2b.S', 'gata3.L',
  'tal1.L', 'nfia.S', 'nfib.L', 'pitx2.S.1'
)

DotPlot(
  vent,
  features = rev(genes_vent),
  assay = 'RNA',
  dot.scale = 9.5
) +
  RotatedAxis() +
  coord_flip() +
  scale_colour_gradient2(low = "white", mid = "gray85", high = "black")

vent <- RenameIdents(
  vent,
  '13' = 'E18', '16' = 'E19', '10' = 'E20', '19' = 'E21', '14' = 'E22',
  '0' = 'E23', '6' = 'E24', '5' = 'E25', '9' = 'E26', '2' = 'M1',
  '7' = 'M2', '3' = 'M3', '11' = 'I1', '17' = 'I2', '12' = 'I3',
  '8' = 'I4', '18' = 'I5', '1' = 'I6', '15' = 'I7', '4' = 'I8',
  '20' = 'I9'
)


## 10. Dorsal inhibitory populations       ##


inh <- subset(scint1, idents = c(4, 1, 2))
inh <- RunPCA(inh, npcs = 50)
inh <- RunUMAP(inh, dims = 1:50)
inh <- FindNeighbors(inh, dims = 1:50)
inh <- FindClusters(inh, resolution = 0.7)

genes_inh_dors <- unique(rev(c(
  'gad2.L', 'slc6a5.S', 'sall3.L', 'adamts5.S',
  'pax2.L', 'nrp1.L', 'ntn1.L', 'nxph4.S', 'lhx1.L',
  'pax5.L', 'LOC108699226', 'LOC108716396', 'tcf12.L',
  'mafb.L', 'pcp4.S', 'nos1.L', 'rspo2.L', 'rorb.L',
  'pyy.L', 'sox5', 'ecel1.L', 'LOC108707572', 'maf.S',
  'asic4.S', 'pex5l.L', 'prox1.L', 'neurod1.S', 'eya4.L'
)))

DotPlot(inh, assay = 'RNA', features = genes_inh_dors) +
  RotatedAxis() +
  coord_flip() +
  scale_colour_gradient2(low = "white", mid = "gray85", high = "black")

levels(inh) <- c(9, 7, 18, 6, 2, 0, 1, 17, 5, 4, 12, 16, 10, 13, 14, 3, 15, 11, 8)

inh <- RenameIdents(
  inh,
  '9' = 'I10', '7' = 'I11', '18' = 'I12', '6' = 'I13', '2' = 'I14',
  '0' = 'I15', '1' = 'I16', '17' = 'I17', '5' = 'I18', '4' = 'I19',
  '12' = 'I20', '16' = 'I21', '10' = 'I22', '13' = 'I23', '14' = 'I24',
  '3' = 'I25', '15' = 'I26', '11' = 'I27', '8' = 'I28'
)

## 11. Combine labels into scint (all neurons data object)         ##


## small cluster 9 = dI2-like excitatory -> E27, and rename MNs
scint <- RenameIdents(scint, '9' = 'E27', '7' = 'MN')

scint$labels <- as.character(Idents(scint))

scint$labels[Cells(exc)]  <- paste(Idents(exc))
scint$labels[Cells(vent)] <- paste(Idents(vent))
scint$labels[Cells(inh)]  <- paste(Idents(inh))

Idents(scint) <- scint$labels

levels(scint) <- c(
  'E1','E2','E3','E4','E5','E6','E7','E8','E9','E10','E11','E12','E13','E14','E15',
  'E16','E17','E18','E19','E20','E21','E22','E23','E24','E25','E26','E27',
  'M1','M2','M3',
  'I1','I2','I3','I4','I5','I6','I7','I8','I9','I10','I11','I12','I13','I14',
  'I15','I16','I17','I18','I19','I20','I21','I22','I23','I24','I25','I26','I27','I28',
  'Pitx2','MN'
)


##clean up metadata and save the object
saveRDS(scint,  '/.../AllCells_Integrated_Xenopus_Adult.rds')


## 12. DotPlot of markers across final NTs ##


genes_all <- rev(unique(c(
  'LOC108713448', 'slc17a7.L', 'gad2.L','slc6a5.S', 'nrxn3.L',
  'slc18a3.S', 'slc10a4.L', 'lmx1b.1.S', 'cpne4.L', 'bnc2.S',
  'nmu.L.1', 'skor1.L','npy.L', 'tac1.L', 'satb1.L', 'cdh23.L',
  'tcf4.S', 'penk.L', 
  'zfhx3.L','zfhx4.L','pou2f2.L', 'prdm16.S','lhx2.S', 'lhx9.S',
  'vsx2.S', 'sim1.S','hmx3.L','foxp2.L', 'foxp1.S', 'esrrg.S',
  'hoxc10.L','hoxc9.L', 'irx5.L','nr4a2.L', 'en1.L', 'tfap2b.S',
  'gata3.L','tal1.L', 'nfia.S', 'nfib.L', 
  'pax2.L', 'ntn1.L', 'nxph4.S','lhx1.L', 'sall3.L', 'tcf12.L', 
  'mafb.L','nos1.L', 'rspo2.L', 'rorb.L', 'sox5', 'ecel1.L',
  'LOC108707572', 'maf.S', 'asic4.S', 'pex5l.L', 'prox1.L', 
  'LOC108696296','pitx2.S.1', 'calca.L', 'chodl.L'
)))

DotPlot(scint, assay = 'RNA', features = genes_all) +
  RotatedAxis() +
  coord_flip() +
  scale_colour_gradient2(low = "white", mid = "gray85", high = "black")


## 13. DimPlot by detailed labels          ##


cols_lab <- c(
  'E1' = 'dodgerblue', 'E2' = 'dodgerblue','E3' = 'dodgerblue','E4' = 'dodgerblue',
  'E5' = 'dodgerblue','E6' = 'dodgerblue','E7' = 'dodgerblue','E8' = 'dodgerblue',
  'E9' = 'dodgerblue','E10' = 'dodgerblue','E11' = 'dodgerblue','E12' = 'dodgerblue',
  'E13' = 'dodgerblue','E14' = 'dodgerblue','E15' = 'dodgerblue','E16' = 'dodgerblue',
  'E17' = 'dodgerblue','E18' = 'sienna1','E19' = 'sienna1','E20' = 'sienna1',
  'E21' = 'sienna1','E22' = 'sienna1','E23' = 'sienna1','E24' = 'sienna1',
  'E25' = 'sienna1','E26' = 'sienna1','E27' = 'sienna1','E28' = 'sienna1',
  'E29' = 'sienna1','E30' = 'sienna1', 
  'I1' = 'sienna1','I2' = 'sienna1','I3' = 'sienna1','I4' = 'sienna1',
  'I5' = 'sienna1','I6' = 'sienna1','I7' = 'sienna1','I8' = 'sienna1',
  'I9' = 'sienna1','I10' = 'sienna1','I11' = 'sienna1','I12' = 'magenta3',
  'I13' = 'magenta3','I14' = 'magenta3','I15' = 'magenta3','I16' = 'magenta3',
  'I17' = 'magenta3','I18' = 'magenta3','I19' = 'magenta3','I20' = 'magenta3',
  'I21' = 'magenta3','I22' = 'magenta3','I23' = 'magenta3','I24' = 'magenta3',
  'I25' = 'magenta3','I26' = 'magenta3','I27' = 'magenta3','I28' = 'magenta3',
  'I29' = 'magenta3','I30' = 'magenta3', 'Pitx2' = 'seagreen', 'MN' = 'green4'
)

DimPlot(scint, label = TRUE, repel = TRUE, label.size = 7, cols = cols_lab) +
  NoLegend() +
  NoAxes()


## 14. Collapse to NT classes and plot     ##



scint <- RenameIdents(
  scint,
  'E1' = 'Exc', 'E2' = 'Exc', 'E3' = 'Exc', 'E4' = 'Exc', 'E5' = 'Exc', 'E6' = 'Exc',
  'E7' = 'Exc', 'E8' = 'Exc', 'E9' = 'Exc', 'E10' = 'Exc', 'E11' = 'Exc',
  'E12' = 'Exc', 'E13' = 'Exc', 'E14' = 'Exc', 'E15' = 'Exc', 'E16' = 'Exc',
  'E17' = 'Exc', 'E18' = 'Exc', 'E19' = 'Exc', 'E20' = 'Exc', 'E21' = 'Exc',
  'E22' = 'Exc', 'E23' = 'Exc', 'E24' = 'Exc', 'E25' = 'Exc', 'E26' = 'Exc',
  'E27' = 'Exc', 'E28' = 'Exc', 'E29' = 'Exc', 'E30' = 'Exc',
  'I1' = 'Inh', 'I2' = 'Inh', 'I3' = 'Inh', 'I4' = 'Inh', 'I5' = 'Inh', 'I6' = 'Inh',
  'I7' = 'Inh', 'I8' = 'Inh', 'I9' = 'Inh', 'I10' = 'Inh', 'I11' = 'Inh',
  'I12' = 'Inh', 'I13' = 'Inh', 'I14' = 'Inh', 'I15' = 'Inh', 'I16' = 'Inh',
  'I17' = 'Inh', 'I18' = 'Inh', 'I19' = 'Inh', 'I20' = 'Inh', 'I21' = 'Inh',
  'I22' = 'Inh', 'I23' = 'Inh', 'I24' = 'Inh', 'I25' = 'Inh', 'I26' = 'Inh',
  'I27' = 'Inh', 'I28' = 'Inh', 'I29' = 'Inh', 'I30' = 'Inh',
  'MN' = 'Cholinergic', 'Pitx2' = 'Cholinergic'
)

pdf('dimplot_scint_by_nts1.pdf', width = 5, height = 5)
colors_nts <- c("Cholinergic" = "lawngreen", "Exc" = "red2", "Inh" = "blue2")
DimPlot(scint, label = FALSE, cols = colors_nts) + NoLegend() + NoAxes()
dev.off()
