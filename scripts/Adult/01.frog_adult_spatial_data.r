library(dplyr)
library(data.table)  
library(tidyr)
library(purrr) 


base_dir <- 'AdultSpatialData'

wrangled_files <- list.files(path = base_dir, pattern = "_wrangled.tsv$", full.names = TRUE, recursive = TRUE)
cell_stats_files <- list.files(path = base_dir, pattern = "_cell_stats.csv$", full.names = TRUE, recursive = TRUE)


#Count matrix
merged_wrangled <- lapply(wrangled_files, function(file) {
  read.csv(file, sep = "\t")
}) %>% Reduce(function(x, y) full_join(x, y, by = "gene"), .)

rownames(merged_wrangled) <- merged_wrangled$gene
merged_wrangled <- merged_wrangled[,-1]

##replace NAs
library(future.apply)

plan(multisession)

spatial <- merged_wrangled
spatial <- as.data.frame(future_lapply(spatial, function(col) {
  col[is.na(col)] <- 0
  col
}))

rownames(spatial) <- rownames(merged_wrangled)

rowname_mapping <- c(
  'chat.L'   = 'LOC108696296',
  'meis2.L'  = 'LOC108699226',
  'npffr2.L' = 'LOC108701554',
  'pou6f2.L' = 'LOC108719010',
  'sox5.L'   = 'sox5',
  'tcf4.L'   = 'Xelaev18011072m',
  'zfpm2.L'  = 'LOC108719377'
)



current_rownames <- rownames(spatial)
new_rownames <- current_rownames
idx <- current_rownames %in% names(rowname_mapping)
new_rownames[idx] <- rowname_mapping[current_rownames[idx]]
rownames(spatial) <- new_rownames


## colnames: replace '.' with '-' to match metadata
colnames(spatial) <- gsub("\\.", "-", colnames(spatial))

## replace NAs with 0
spatial[is.na(spatial)] <- 0


## Spatial object and coarse LT       

spatial <- CreateSeuratObject(spatia)
spatial <- SCTransform(spatial)

##load all cells single-nuclei adult frog object

all <- readRDS('/.../AllCells_Integrated_Xenopus_Adult.rds')

DefaultAssay(all) <- 'RNA'
all <- JoinLayers(all)
all <- NormalizeData(all)
all <- ScaleData(all, features = rownames(spatial1))
all <- RunPCA(all, features = rownames(spatial1))
all <- RunUMAP(all, dims = 1:10)

## SCTransform with spatial genes as residual features
all <- SCTransform(all, residual.features = rownames(spatial1))
all <- RunPCA(all)
all <- RunUMAP(all, dims = 1:5)


## collapse to Neurons / others
all <- RenameIdents(all, 'INs' = 'Neurons', 'MNs' = 'Neurons')

all <- JoinLayers(all)
all <- GetAssayData(all, assay = 'RNA', slot = 'counts')
all <- SCTransform(all, variable.features.n = 10000)


##Label transfer to get neurons in spatial data               

anchors <- FindTransferAnchors(
  reference = all,
  query = spatial,
  normalization.method = "SCT"
)

predictions.assay <- TransferData(
  anchorset = anchors,
  refdata = Idents(all),
  prediction.assay = TRUE
)

spatial[['predictions']] <- predictions.assay

spatial$prediction.id <- GetTransferPredictions(spatial, score.filter = 0.5)
Idents(spatial) <- spatial$prediction.id

##clean up meta data and save spatial all cells object
saveRDS(spatial, '/.../Xenopus_Adult_SpatialData_AllCells.rds')


#subset neurons
neu_sp <- subset(spatial, idents = 'Neurons')


## neurotransmitter-profile (exc, inh, mns) - level label transfer to spatial ##              

scint <- subset(scint, idents = c(levels(scint1)[1:60], 'MN'))
scint$labels <- Idents(scint)

scint <- RenameIdents(
  scint,
  'E1' = 'Exc', 'E2' = 'Exc', 'E3' = 'Exc', 'E4' = 'Exc', 'E5' = 'Exc',
  'E6' = 'Exc', 'E7' = 'Exc', 'E8' = 'Exc', 'E9' = 'Exc', 'E10' = 'Exc', 
  'E11' = 'Exc', 'E12' = 'Exc', 'E13' = 'Exc', 'E14' = 'Exc', 'E15' = 'Exc',
  'E16' = 'Exc', 'E17' = 'Exc', 'E18' = 'Exc', 'E19' = 'Exc', 'E20' = 'Exc',
  'E21' = 'Exc', 'E22' = 'Exc', 'E23' = 'Exc', 'E24' = 'Exc', 'E25' = 'Exc',
  'E26' = 'Exc', 'E27' = 'Exc', 'E28' = 'Exc', 'E29' = 'Exc', 'E30' = 'Exc',
  'I1' = 'Inh', 'I2' = 'Inh', 'I3' = 'Inh', 'I4' = 'Inh', 'I5' = 'Inh',
  'I6' = 'Inh', 'I7' = 'Inh', 'I8' = 'Inh', 'I9' = 'Inh', 'I10' = 'Inh',
  'I11' = 'Inh', 'I12' = 'Inh', 'I13' = 'Inh', 'I14' = 'Inh', 'I15' = 'Inh',
  'I16' = 'Inh', 'I17' = 'Inh', 'I18' = 'Inh', 'I19' = 'Inh', 'I20' = 'Inh',
  'I21' = 'Inh', 'I22' = 'Inh', 'I23' = 'Inh', 'I24' = 'Inh', 'I25' = 'Inh',
  'I26' = 'Inh', 'I27' = 'Inh', 'I28' = 'Inh', 'I29' = 'Inh', 'I30' = 'Inh'
)

scint$nts <- as.character(Idents(scint))

## prep neurons_sn for LT
neurons_sn <- scint
DefaultAssay(neurons_sn) <- 'RNA'
neurons_sn <- JoinLayers(neurons_sn)
neurons_sn <- GetAssayData(neurons_sn, assay = 'RNA', slot = 'counts')
neurons_sn <- CreateSeuratObject(neurons_sn)
neurons_sn <- SCTransform(neurons_sn, variable.features.n = 10000)
Idents(neurons_sn) <- Idents(scint)

neu_sp <- SCTransform(neu_sp)

anchors <- FindTransferAnchors(
  reference = neurons_sn,
  query = neu_sp,
  normalization.method = "SCT"
)

predictions.assay <- TransferData(
  anchorset = anchors,
  refdata = Idents(neurons_sn),
  prediction.assay = TRUE
)

neu_sp[['predictions']] <- predictions.assay
neu_sp$prediction.id <- GetTransferPredictions(neu_sp, score.filter = 0.5)
Idents(neu_sp) <- neu_sp$prediction.id
neu_sp$nts <- Idents(neu_sp)


##  Excitatory families on spatial


Idents(scint) <- scint$labels
scint_exc <- subset(scint, idents = levels(scint)[1:30])
DefaultAssay(scint_exc) <- 'RNA'
scint_exc <- JoinLayers(scint_exc)
exc_labels <- Idents(scint_exc)
scint_exc <- GetAssayData(scint_exc, assay = 'RNA', slot = 'counts')
scint_exc <- CreateSeuratObject(scint_exc)
scint_exc <- SCTransform(scint_exc, variable.features.n = 10000)
Idents(scint_exc) <- exc_labels


## SCTransform with spatial genes as residual features for grouping
scint_exc <- SCTransform(scint_exc, residual.features = rownames(neu_sp))
scint_exc <- RunPCA(scint_exc)
scint_exc <- RunUMAP(scint_exc, dims = 1:10)
scint_exc$labels <- Idents(scint_exc)
scint_exc <- FindNeighbors(scint_exc, dims = 1:10)
scint_exc <- FindClusters(scint_exc, resolution = 2)


scint_exc$labels <- Idents(scint_exc)
scint_exc <- RenameIdents(
  scint_exc,
  'E1' = 'Exc_Cpne4/Bnc2', 'E2' = 'Exc_Cpne4/Bnc2', 'E3' = 'Exc_Cpne4/Bnc2',
  'E4' = 'Exc_Cpne4/Bnc2', 'E5' = 'Exc_Cpne4/Bnc2', 'E6' = 'Exc_Cpne4/Bnc2',
  'E7' = 'Exc_Cpne4/Bnc2', 'E8' = 'Exc_Skor1', 'E9' = 'Exc_Skor1',
  'E10' = 'Exc_Skor1', 'E11' = 'Exc_Npy', 'E12' = 'Exc_Npy',
  'E13' = 'Exc_Tcf4/Penk', 'E14' = 'Exc_Tcf4/Penk', 'E15' = 'Exc_Tcf4/Penk',
  'E16' = 'Exc_Tcf4/Penk', 'E17' = 'Exc_Tcf4/Penk',
  'E18' = 'Exc_Vent', 'E19' = 'Exc_Vent', 'E20' = 'Exc_Vent', 'E21' = 'Exc_Vent',
  'E22' = 'Exc_Vent', 'E23' = 'Exc_Vent', 'E24' = 'Exc_Vent', 'E25' = 'Exc_Vent',
  'E26' = 'Exc_Vent', 'E27' = 'Exc_Vent', 'E28' = 'Exc_Vent',
  'E29' = 'Exc_Vent', 'E30' = 'Exc_Vent'
)

exc_sp <- subset(neu_sp, idents = 'Exc')
exc_sp <- SCTransform(exc_sp)

anchors <- FindTransferAnchors(
  reference = scint_exc,
  query = exc_sp,
  normalization.method = "SCT"
)

predictions.assay <- TransferData(
  anchorset = anchors,
  refdata = Idents(scint_exc),
  prediction.assay = TRUE
)

exc_sp[['predictions']] <- predictions.assay
exc_sp$prediction.id <- GetTransferPredictions(exc_sp, score.filter = 0.3)
Idents(exc_sp) <- exc_sp$prediction.id


## spatial plots and densities for excitatory families
Idents(exc_sp) <- exc_sp$slide
c2_3_exc_sp <- subset(exc_sp, idents = 'C2_3')
c2_3_exc_sp <- subset(c2_3_exc_sp, subset = y < 10000)
c2_3_exc_sp <- subset(c2_3_exc_sp, subset = x > 1000)
Idents(c2_3_exc_sp) <- c2_3_exc_sp$prediction.id
c2_3_exc_sp <- subset(c2_3_exc_sp, idents = levels(c2_3_exc_sp)[1:5])

ggplot(c2_3_exc_sp@meta.data, aes(x = x, y = y, color = prediction.id)) +
  geom_point(alpha = 0.8, stroke = 0, size = 2.4) +
  scale_color_manual(values = c(
    "Exc-Vent" = "tomato",
    "Exc-Cpne4/Bnc2" = 'lightskyblue',
    "Exc-Npy" = 'cornflowerblue',
    "Exc-Tcf4/Penk" = "navy",
    "Exc-Skor1" = 'blue'
  )) +
  theme_void() +
  theme(legend.position = "none")


meta_data <- c2_3_exc_sp@meta.data
subset_data <- subset(meta_data, y < 10000)
subset_data <- subset(subset_data, x > 1000)


##density plots
density_data <- subset_data %>%
  group_by(prediction.id) %>%
  do(data.frame(
    density = density(.$x)$y,
    x       = density(.$x)$x
  ))

density_data <- density_data %>%
  group_by(prediction.id) %>%
  mutate(normalized_density = density / max(density))

custom_colors <- c(
  "Exc-Vent" = "tomato",
  "Exc-Cpne4/Bnc2" = 'lightskyblue',
  "Exc-Npy" = 'cornflowerblue',
  "Exc-Tcf4/Penk" = "navy",
  "Exc-Skor1" = 'blue'
)

              
ggplot(density_data, aes(x = x, y = normalized_density,
                         color = prediction.id, fill = prediction.id)) +
  geom_ribbon(aes(ymin = 0, ymax = normalized_density), alpha = 0.3) +
  geom_line() +
  labs(
    title = "Normalized Density of Cell Distributions by Prediction ID",
    x = "X-axis",
    y = "Normalized Density"
  ) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  theme(legend.position = "right")


##  Inhibitory families on spatial      ##
scint_inh <- subset(scint, idents = levels(scint)[31:60])
DefaultAssay(scint_inh) <- 'RNA'
scint_inh <- JoinLayers(scint_inh)
inh_labels <- Idents(scint_inh)
scint_inh <- GetAssayData(scint_inh, assay = 'RNA', slot = 'counts')
scint_inh <- CreateSeuratObject(scint_inh)
scint_inh <- SCTransform(scint_inh, variable.features.n = 10000)
Idents(scint_inh) <- inh_labels


scint_inh <- SCTransform(scint_inh, residual.features = rownames(neu_sp))
scint_inh <- RunPCA(scint_inh)
scint_inh <- RunUMAP(scint_inh, dims = 1:10)
scint_inh$labels <- Idents(scint_inh)
scint_inh <- FindNeighbors(scint_inh, dims = 1:10)
scint_inh <- FindClusters(scint_inh, resolution = 2)

scint_inh$labels <- Idents(scint_inh)
scint_inh <- RenameIdents(
  scint_inh,
  'I1' = 'Inh_Vent', 'I2' = 'Inh_Vent', 'I3' = 'Inh_Vent', 'I4' = 'Inh_Vent',
  'I5' = 'Inh_Vent', 'I6' = 'Inh_Vent', 'I7' = 'Inh_Vent', 'I8' = 'Inh_Vent',
  'I9' = 'Inh_Vent', 'I10' = 'Inh_Vent', 'I11' = 'Inh_Vent',
  'I12' = 'Inh_Ntn1', 'I13' = 'Inh_Ntn1',
  'I14' = 'Inh_Other', 'I15' = 'Inh_Other', 'I16' = 'Inh_Other',
  'I17' = 'Inh_Tcf4/Nfi', 'I18' = 'Inh_Tcf4/Nfi',
  'I19' = 'Inh_Nos1', 'I20' = 'Inh_Nos1',
  'I21' = 'Inh_Other', 'I22' = 'Inh_Other',
  'I23' = 'Inh_Rorb', 'I24' = 'Inh_Rorb',
  'I25' = 'Inh_Tcf4/Nfi', 'I26' = 'Inh_Tcf4/Nfi',
  'I27' = 'Inh_Tcf4/Nfi', 'I28' = 'Inh_Tcf4/Nfi',
  'I29' = 'Inh_Tcf4/Nfi', 'I30' = 'Inh_Tcf4/Nfi'
)

inh_sp <- subset(neu_sp, idents = 'Inh')
inh_sp <- SCTransform(inh_sp)

anchors <- FindTransferAnchors(
  reference = scint_inh,
  query = inh_sp,
  normalization.method = "SCT"
)

predictions.assay <- TransferData(
  anchorset = anchors,
  refdata = Idents(scint_inh),
  prediction.assay = TRUE
)

inh_sp[['predictions']] <- predictions.assay
inh_sp$prediction.id <- GetTransferPredictions(inh_sp, score.filter = 0.3)
Idents(inh_sp) <- inh_sp$prediction.id


Idents(inh_sp) <- inh_sp$slide
c2_3_inh_sp <- subset(inh_sp, idents = 'C2_3')
c2_3_inh_sp <- subset(c2_3_inh_sp, subset = y < 10000)
c2_3_inh_sp <- subset(c2_3_inh_sp, subset = x > 1000)
Idents(c2_3_inh_sp) <- c2_3_inh_sp$prediction.id
c2_3_inh_sp <- subset(c2_3_inh_sp, idents = levels(c2_3_inh_sp)[1:6])


ggplot(c2_3_inh_sp@meta.data, aes(x = x, y = y, color = prediction.id)) +
  geom_point(alpha = 0.8, stroke = 0, size = 3.2) +
  scale_color_manual(values = c(
    "Inh-Vent" = "orange", "Inh-Nos1" = 'hotpink',
    "Inh-Ntn1" = 'purple', "Inh-Other" = "maroon",
    "Inh-Tcf4/Nfi" = 'magenta4', "Inh-Rorb" = 'violet'
  )) +
  theme_void() +
  theme(legend.position = "none")


meta_data <- c2_3_inh_sp@meta.data
subset_data <- subset(meta_data, y < 10000)
subset_data <- subset(subset_data, x > 1000)

density_data <- subset_data %>%
  group_by(prediction.id) %>%
  do(data.frame(
    density = density(.$x)$y,
    x       = density(.$x)$x
  )) %>%
  group_by(prediction.id) %>%
  mutate(normalized_density = density / max(density))

custom_colors <- c(
  "Inh-Vent" = "orange", "Inh-Nos1" = 'hotpink',
  "Inh-Ntn1" = 'purple', "Inh-Other" = "maroon",
  "Inh-Tcf4/Nfi" = 'magenta4', "Inh-Rorb" = 'violet'
)


ggplot(density_data, aes(x = x, y = normalized_density,
                         color = prediction.id, fill = prediction.id)) +
  geom_ribbon(aes(ymin = 0, ymax = normalized_density), alpha = 0.3) +
  geom_line() +
  labs(
    title = "Normalized Density of Cell Distributions by Prediction ID",
    x = "X-axis",
    y = "Normalized Density"
  ) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  theme(legend.position = "right")

# Final combined spatial neu object   ##


Idents(neu_sp) <- neu_sp$prediction.id
neu_sp$groups <- as.character(Idents(neu_sp))
neu_sp$groups[Cells(exc_sp)] <- paste(exc_sp$prediction.id)
neu_sp$groups[Cells(inh_sp)] <- paste(inh_sp$prediction.id)

Idents(neu_sp) <- neu_sp$groups
neu_sp <- subset(neu_sp, idents = c(levels(neu_sp)[1:9], levels(neu_sp)[11:13]))

levels(neu_sp) <- c(
  "Exc-Skor1","Exc-Cpne4/Bnc2", "Exc-Npy", "Exc-Tcf4/Penk",
  "Exc-Vent","Inh-Vent", "Inh-Ntn1","Inh-Other", "Inh-Nos1",
  "Inh-Rorb", "Inh-Tcf4/Nfi", "MN"
)

neu_sp$groups <- Idents(neu_sp)


DotPlot(
  neu_sp,
  features = rev(unique(c(
    'slc17a7.L','slc6a5.S','slc18a3.S','lmx1b.1.S','cpne4.L',
    'bnc2.S','prkcg.S',"kcnq4.L",'rasal1.S','npy.L','tac1.L',
    'cdh23.L','penk.L','Xelaev18011072m','LOC108699226','zfhx4.L',
    'zfhx3.L','foxp2.L','LOC108719377','ntn1.L','sst.2','sall3.L',
    'nos1.L','mafb.L','rorb.L','maf.L','nfib.L','neurod2.S',
    "adamts5.S",'ntn1.L',"slit2.S"
  ))),
  dot.scale = 8,
  assay = 'RNA',
  cols = c('white', 'black')
) +
  RotatedAxis() +
  coord_flip()
dev.off()

##cleaning up metadata and saving the object
saveRDS(neu_sp, '/.../Xenopus_Adult_SpatialData_Neurons.rds')

