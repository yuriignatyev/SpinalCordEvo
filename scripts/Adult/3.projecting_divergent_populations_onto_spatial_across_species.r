# Here, we are determining the spatial location of divergent vs conserved populations defined in cross-species analyses
# We first rename clusters in the integrated subsetted data of dorsal excitatory and inhibitory neurons into 'Divergent' and 'Conserved' based on species mixing (see Fig.4G and S18)
# And then use these assignments to project them into frog spatial data and mouse Visium data from Kathe et al.

library(Seurat)
library(ggplot2)

#Let's first load all the objects we will need

#Individual neuronal snRNA-seq data objects - we will need frog and mouse (Russ et al.) data
#mouse
sc.list <- readRDS('/.../List_of_Individual_SpeciesStudy_Neuronal_Objects_FrogMouseHumanZebrafish_SCTrasformed_Ready_for_Integration.rds')
russ <- sc.list[[2]]
#frog 
frog <- readRDS('/.../Neurons_Xenopus_Adult.rds')

#Integrated data 
fmh <- readRDS('/.../Neurons_Integrated_Across_FrogMouseHuman_CCA.rds')
exc <- readRDS('/.../DorsalExcitatoryNeurons_Integrated_Subset_Across_FrogMouseHuman_CCA.rds')
inh <- readRDS('/.../DorsalInhibitoryNeurons_Integrated_Subset_Across_FrogMouseHuman_CCA.rds')

#Spatial datasets (here we'll use already processed neuronal objects with projections in 'predictions' assay, but we'll re-project)
frog_sp <- readRDS('/.../Xenopus_Adult_SpatialData_Neurons.rds')
visium <- readRDS('/.../Kathe2022_MouseVisiumSpatialData_Neurons.rds')

#1. Projecting onto the frog spatial dataset and visualizing

#Prepping frog snRNA-seq data for label transfer
DefaultAssay(frog) <- 'RNA'
frog <- JoinLayers(frog)
frog <- GetAssayData(frog, assay = 'RNA', slot = 'counts')
frog <- CreateSeuratObject(frog)
frog <- SCTransform(frog, variable.features.n = 10000)

#Relabelling clusters in integrated data into 'Divergent' and 'Conserved' and taking these labels into the original frog data ('frog')
Idents(exc) <- exc$species
exc_frog <- subset(exc, idents = 'Frog')
Idents(exc_frog) <- exc_frog$seurat_clusters
exc_frog <- RenameIdents(exc_frog, '11' = 'Divergent', '12' = 'Divergent', '0' = 'Conserved', '1' = 'Conserved','2' = 'Conserved','3' = 'Conserved','4' = 'Conserved','5' = 'Conserved','6' = 'Conserved','7' = 'Conserved','8' = 'Conserved','9' = 'Conserved','10' = 'Conserved','13' = 'Conserved','14' = 'Conserved','15' = 'Conserved','16' = 'Conserved','18' = 'Conserved','19' = 'Conserved','20' = 'Conserved','22' = 'Conserved')

Idents(fmh) <- fmh$species
frog_int <- subset(fmh, idents = 'Frog')

frog_int$consdiv <- 'Conserved'
frog_int$lab[Cells(exc_frog)] <- paste(Idents(exc_frog))
frog@meta.data$consdiv <- frog_int@meta.data$consdiv

Idents(frog) <- frog$consdiv

#Now that we have 'Divergent' and 'Conserved' labels inthe  original frog snRNA-seq data, let's project it onto the spatial frog dataset via label transfer
anchors <- FindTransferAnchors(reference = frog, query = frog_sp, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = Idents(frog), prediction.assay = TRUE)

frog_sp[['predictions']]<- predictions.assay

#Let's visualize on frog spatial data

DefaultAssay(frog_sp) <- "predictions"

Idents(frog_sp) <- frog_sp$slide

c2_3_frog_sp <- subset(frog_sp, idents = c('C2_3'))
c2_3_frog_sp <- subset(c2_3_frog_sp, subset = y < 10000)
c2_3_frog_sp <- subset(c2_3_frog_sp, subset = x > 1000)

df <- FetchData(c2_3_frog_sp, vars = c("x","y","Conserved","Divergent"))
df_plot <- df[order(df$Divergent), ]

ggplot(df_plot, aes(x = x, y = y, color = Divergent)) +
  geom_point(size = 1) +
  coord_equal() +
  scale_color_gradient(low = "gray90", high = "black") +
  theme_void() +
  ggtitle("prediction score")


#2. Projecting onto the mouse spatial dataset and visualizing

Idents(exc) <- exc$study  
Idents(inh) <- inh$study  

exc_russ <- subset(exc, idents = 'Russ')
inh_russ <- subset(inh, idents = 'Russ')

Idents(exc_russ) <- exc_russ$seurat_clusters
Idents(inh_russ) <- inh_russ$seurat_clusters


exc_russ <- RenameIdents(exc_russ, '0' = 'Divergent','1' = 'Divergent','2' = 'Conserved','3' = 'Divergent','4' = 'Divergent','5' = 'Conserved','6' = 'Conserved','7' = 'Divergent','8' = 'Conserved','9' = 'Divergent','10' = 'Divergent','11' = 'Conserved','12' = 'Conserved','13' = 'Divergent','14' = 'Conserved','15' = 'Conserved','16' = 'Divergent','17' = 'Divergent','18' = 'Divergent','19' = 'Divergent','20' = 'Conserved','21' = 'Divergent', '22' = 'Conserved')
inh_russ <- RenameIdents(inh_russ, '0' = 'Divergent', '1' = 'Conserved', '2' = 'Conserved', '3' = 'Conserved', '4' = 'Conserved', '5' = 'Conserved', '6' = 'Conserved', '7' = 'Conserved', '8' = 'Conserved', '9' = 'Conserved', '10' = 'Conserved', '11' = 'Conserved', '12' = 'Conserved', '13' = 'Conserved', '14' = 'Conserved', '15' = 'Conserved', '16' = 'Conserved', '17' = 'Divergent', '18' = 'Divergent', '19' = 'Divergent', '20' = 'Divergent')

Idents(fmh) <- fmh$study
russ_int <- subset(fmh, idents = 'Russ')

russ_int$consdiv <- 'Conserved'
russ_int$consdiv[Cells(exc_russ)] <- paste(Idents(exc_russ))
russ_int$consdiv[Cells(inh_russ)] <- paste(Idents(inh_russ))

#Taking them into Russ et al. neuronal object
russ@meta.data$consdiv <- as.data.frame(russ_int@meta.data$consdiv)
russ$consdiv <- russ_int@meta.data$consdiv
Idents(russ) <- russ$consdiv

#For mouse, we'll use 'RNA' assay for label transfer since it worked better for label transfer between Russ et al. snRNA_seq data and visium spatial data from Kathe et al.

DefaultAssay(visium) <- 'RNA'
DefaultAssay(russ) <- 'RNA'

russ <- NormalizeData(russ)
visium <- NormalizeData(visium)

russ <- FindVariableFeatures(russ, selection.method = 'vst', nfeatures = 10000)
visium <- FindVariableFeatures(visium, selection.method = 'vst', nfeatures = 10000)

anchors <- FindTransferAnchors(reference = russ, query = visium, normalization.method = "LogNormalize", features = intersect(VariableFeatures(russ), VariableFeatures(visium)))
predictions.assay <- TransferData(anchorset = anchors, refdata = Idents(russ), prediction.assay = TRUE)
visium[['predictions']] <- predictions.assay

#Let's visualize on mouse spatial data

df <- FetchData(visium, vars = c("coord_x","coord_y","Conserved","Divergent"))
df_plot <- df[order(df$Divergent), ]

ggplot(df_plot, aes(x = coord_x, y = coord_y, color = Divergent)) +
  geom_point(size = 1) +
  coord_equal() +
  scale_color_gradient(low = "gray95", high = "black") +
  theme_void() +
  ggtitle("prediction score")


