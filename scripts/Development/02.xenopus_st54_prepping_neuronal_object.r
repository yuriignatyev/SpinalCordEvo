###subsetting interneurons and motoneurons from 'all cells' object 'AllCells_Filtered_Labelled_Xenopus54_After_MixtoolsTresholding_and_RemovingDoublets.rds'
library(Seurat)
library(ggplot2)
library(dplyr)

## ---------- load ----------
all <- readRDS('/.../DevelopmentSingleCell/AllCells_Filtered_Labelled_Xenopus54_After_MixtoolsTresholding_and_RemovingDoublets.rds')

levels(all)

## ---------- subset neurons ----------
frog_neu <- subset(all, idents = levels(all)[5:13])

frog_neu <- NormalizeData(frog_neu, verbose = FALSE)
frog_neu <- FindVariableFeatures(frog_neu, selection.method = 'vst', nfeatures = 10000)
frog_neu <- ScaleData(frog_neu)
frog_neu <- RunPCA(frog_neu, npcs = 50)
frog_neu <- RunUMAP(frog_neu, dims = 2:50)

frog_neu <- FindNeighbors(frog_neu, dims = 2:50, verbose = FALSE)

##clusterting frog and mouse gives almost the same number of clusters 
### we used high resolutions to get more precise labels since some of the neural classes can have overlapping expression
#-wise but still expose important biological differences on a levels of a small number genes

frog_neu <- FindClusters(frog_neu, resolution = 2)

##let's see how cardinal class markers map onto these clusters
##we also did differential expression analysis on each cluster and confirmed cardinal class markers being top in diff.exp. markers in most clusters
DotPlot(frog_neu, assay = 'RNA', cols = c('white', 'black'), features = c('btg2.S',
                                                          'nhlh1.L', 'barhl1.L', 'barhl2.L', 'lhx2.S', 'lhx9.S',
                                                          'hmx3.L', 'foxd3.L', 'isl1.L.1', 'tlx3.L', 'drgx.L', 'lmx1b.1.S',
                                                          'lbx1.L', 'sall3.L', 'pax2.L', 'pax8.L', 'wt1.L', 'dmrt3.L', 'evx1.L',
                                                          'evx2.S', 'pitx2.S.1', 'en1.L', 'foxp2.L', 'lhx3.L', 'vsx2.S', 'sox14.L', 
                                                          'sox21.L', 'gata2.L', 'gata3.L', 'tal1.L', 'nkx2-2.L', 'sim1.S', 'slc18a3.S', 
                                                          'slit2.S'
                                                         )) + 
                                                                  RotatedAxis() + coord_flip()


###cluster 13 has low quality cells, so lets remove it and recluster
thirteen <- subset(frog_neu, idents = '13')
thirteen <- RunPCA(thirteen, npcs = 50)
thirteen <- RunUMAP(thirteen, dims = 1:30)
DimPlot(thirteen)
thirteen <- FindNeighbors(thirteen, dims = 1:30)
thirteen <- FindClusters(thirteen, resolution = 0.3)
DimPlot(thirteen, label = TRUE)
VlnPlot(thirteen, c('nFeature_RNA'))
thirteen <- RenameIdents(thirteen, '0' = 'low', '1' = 'good', '2' = 'good')
markers_bad_cells <- FindMarkers(thirteen, ident.1 = 'low', min.pct = 0.25)
FeaturePlot(thirteen, order = TRUE, cols = c('white', 'black'), features = c('nFeature_RNA'))


frog_neu$removing <- as.character(Idents(frog_neu))
frog_neu$removing[Cells(thirteen)] <- paste(Idents(thirteen))
Idents(frog_neu) <- frog_neu$removing
clusters_to_keep <- c(levels(frog_neu)[1:32], levels(frog_neu)[34:47]) 

##lets remove last low quality cells from neuronal object
frog_neu <- subset(frog_neu, idents = clusters_to_keep)
frog_neu <- RunPCA(frog_neu, npcs = 50)
frog_neu <- RunUMAP(frog_neu, dims = 1:50)
frog_neu <- FindNeighbors(frog_neu, dims = 1:50)
frog_neu <- FindClusters(frog_neu, resolution = 2)

## '0' cluster consisted mostly of dI1 cells but also some lower quality cells so lets remove it with reclustering

zero <- subset(frog_neu, idents = 0)
zero <- RunPCA(zero, npcs = 50)
zero <- RunUMAP(zero, dims = 1:10)
zero <- FindNeighbors(zero, dims = 1:10)
zero <- FindClusters(zero, resolution = 0.3)
VlnPlot(zero, c('nFeature_RNA')) + FeaturePlot(zero, c('nFeature_RNA'), cols = c('white', 'black')) + DimPlot(zero, label = TRUE)
zero <- RenameIdents(zero, '0' = 'bad', '1' = 'good', '2' = 'good', '3' = 'good', '4' = 'good')

frog_neu$removing2 <- as.character(Idents(frog_neu))
frog_neu$removing2[Cells(zero)] <- paste(Idents(zero))

##lets remove last low quality cells from neuronal object again
Idents(frog_neu) <- frog_neu$removing2
frog_neu <- subset(frog_neu, idents = c(levels(frog_neu)[1:45], c(46,26)))
frog_neu <- RunPCA(frog_neu, npcs = 50)
frog_neu <- RunUMAP(frog_neu, dims = 2:50)

frog_neu <- FindNeighbors(frog_neu, dims = 2:50)
frog_neu <- FindClusters(frog_neu, resolution = 2)

##cluster 5 has mixed NT identity,so lets take a look
five_markers <- FindMarkers(frog_neu, ident.1 = 5, min.pct = 0.25)

##it is a mix of dI2 and V1 and was marked by foxp4, foxp2, foxd3 which mark multiple classes
five <- subset(frog_neu, idents = 5)
five <- RunPCA(five, npcs = 10)
five <- RunUMAP(five, dims = 1:10)
five <- FindNeighbors(five, dims = 1:10)
five <- FindClusters(five, resolution = 0.2)

DimPlot(five, label = TRUE) + DotPlot(five, features = c('en1.L', 'foxd3.L', 'gad2.L', 'slc17a7.S'))

five <- RenameIdents(five, '0' = 'V1', '1' = 'dI2', '2' = 'V1')

##putting these labels into the neural object

frog_neu$relabelling1 <- as.character(Idents(frog_neu))
frog_neu$relabelling1[Cells(five)] <- paste(Idents(five))

Idents(frog_neu) <- frog_neu$relabelling1

##also cluster 22 had mixed identiteties of dI2, dI1 and V0 excitatory neurons
twentytwo <- subset(frog_neu, idents = 22)
twentytwo <- RunPCA(twentytwo, npcs = 50)
twentytwo <- RunUMAP(twentytwo, dims = 3:10) 
twentytwo <- FindNeighbors(twentytwo, dims = 3:10)
twentytwo <- FindClusters(twentytwo, resolution = 0.2)

DimPlot(twentytwo, label = TRUE) + DotPlot(twentytwo, features = c('evx1.L', 'evx2.S', 'lhx2.S', 'lhx9.S', 'foxd3.L', 'hmx3.L','gad2.L', 'slc17a7.S'))

twentytwo <- RenameIdents(twentytwo, '0' = 'dI2', '1' = 'dI1', '2' = 'V0')

##putting these ones too
frog_neu$relabelling2 <- as.character(Idents(frog_neu))
frog_neu$relabelling2[Cells(twentytwo)] <- paste(Idents(twentytwo))
Idents(frog_neu) <- frog_neu$relabelling2


DotPlot(frog_neu, assay = 'RNA', cols = c('white', 'black'), features = c('btg2.S',
                                                          'nhlh1.L', 'barhl1.L', 'barhl2.L', 'lhx2.S', 'lhx9.S',
                                                          'hmx3.L', 'foxd3.L', 'isl1.L.1', 'tlx3.L', 'drgx.L', 'lmx1b.1.S',
                                                          'lbx1.L', 'sall3.L', 'pax2.L', 'pax8.L', 'wt1.L', 'dmrt3.L', 'evx1.L',
                                                          'evx2.S', 'pitx2.S.1', 'en1.L', 'foxp2.L', 'lhx3.L', 'vsx2.S', 'sox14.L', 
                                                          'sox21.L', 'gata2.L', 'gata3.L', 'tal1.L', 'nkx2-2.L', 'sim1.S', 'slc18a3.S', 
                                                          'slit2.S', 'zfhx3.L', 'nfib.L', 'gad2.L', 'slc17a7.S', 
                                                          'plp1.L'
                                                         )) + 
                                                                     RotatedAxis() + coord_flip()

##finally renaming
frog_neu <- RenameIdents(frog_neu, '0' = 'Differentiating','1' = 'dI4','2' = 'dI5','3' = 'dI4',
'4' = 'dI6','6' = 'dI4','7' = 'dI4','8' = 'dI5','9' = 'V2a','10' = 'V0','11' = 'dI4','12' = 'dI1',
'13' = 'Differentiating','14' = 'V2b','15' = 'MNs','16' = 'dI3','17' = 'V1','18' = 'dI2','19' = 'dI4','20' = 'dI5','21' = 'dI4',
'23' = 'MNs','24' = 'dI5','25' = 'dI4','26' = 'V1','27' = 'dI1','28' = 'dI5',
'29' = 'dI5','30' = 'dI4','31' = 'V2a','32' = 'dI4','33' = 'dI5','34' = 'V1','35' = 'dI5',
'36' = 'V3','37' = '37','38' = 'dI5','39' = 'dI2','40' = 'V3','41' = '41','42' = 'dI4',
'43' = 'Differentiating','44' = 'MNs','45' = 'V0','46' = 'Contamination_Oligodendrocytes')

###clusters 37 and 41 
markers_37 <- FindMarkers(frog_neu, ident.1 = '37', min.pct = 0.25, logfc.threshold = 0.25)
markers_41 <- FindMarkers(frog_neu, ident.1 = '41', min.pct = 0.25, logfc.threshold = 0.25)

##41 are sensory neurons contamination
##37 - no clear cardinal class identity and a lot of Irx genes etc. that are related to Dev processes


cardinal_classes <- c('dI1','dI2','dI3', 'dI4','dI5','dI6','V0','V1','V2a','V2b','V3','MNs')
###subset cardinal classes 
frog_neu <- subset(frog_neu, idents = cardinal_classes)
frog_neu <- RunPCA(frog_neu, npcs = 50)
frog_neu <- RunUMAP(frog_neu, dims = 2:50)

levels(frog_neu) <- cardinal_classes


## ---------- clean up metadata and save final object ----------
saveRDS(frog_neu, '/.../DevelopmentSingleCell/Xenopus54_Neurons_CardinalClasses_Final_Labelled.rds')