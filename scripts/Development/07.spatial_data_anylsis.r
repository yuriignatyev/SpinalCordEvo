library(Seurat)
##Preparing single-cell references for LT onto spatial data

##All cell types
all <- readRDS('/.../AllCells_Filtered_Labelled_Xenopus54_After_MixtoolsTresholding_and_RemovingDoublets.rds')
##Mature Neurons and MNs
neu <- readRDS('/.../Xenopus54_Neurons_CardinalClasses_Final_Labelled.rds')

all <- RenameIdents(all,'Neural Progenitors 2' = 'Neural Progenitors','Neural Progenitors 3' = 'Neural Progenitors',
'Differentiating Neurons 1' = 'Differentiating Neurons','Differentiating Neurons 2' = 'Differentiating Neurons',
'Motor Neurons' = 'Neurons', 'Neural Crest' = 'Other','OCPs' = 'Other','Schwann Cells' = 'Other',
'Mesoderm-derived' = 'Other','Myocytes' = 'Other', 'Immune' = 'Other','Endothelian' = 'Other','Skin' = 'Other')


all$labels <- as.character(Idents(all))
all$labels[Cells(neu)] <- paste(Idents(neu))

Idents(all) <- all$labels

all <- subset(all, idents = c(
  "Interneurons 1", "Peripheral Sensory Neurons", "Interneurons 5",
  "Interneurons 4", "Interneurons 2", "Cerebrospinal Fluid-Contacting Neurons",
  "Serotonergic/Dopaminergic Neurons", "Interneurons 3", "Neurons", "Interneurons 6"
), invert = TRUE)


#separating FP from progenitors
fp <- subset(all, idents = "FP and Neural Progenitors 1")
fp <- FindClusters(fp, resolution = 0.5)
fp <- RenameIdents(fp, '0' = 'Neural Progenitors', '1' = 'Neural Progenitors',
'2' = 'Neural Progenitors','3' = 'Neural Progenitors','4' = 'Neural Progenitors',
'5' = 'Neural Progenitors','6' = 'FP','7' = 'Neural Progenitors')

all$labels <- as.character(Idents(all))
all$labels[Cells(fp)] <- paste(Idents(fp))
Idents(all) <- all$labels


all <- RenameIdents(all, 'MNs' = 'Neurons', 'dI4' = 'Neurons','V1' = 'Neurons',
'dI1' = 'Neurons','dI5' = 'Neurons','V2a' = 'Neurons','dI3' = 'Neurons',
'V0' = 'Neurons','dI2' = 'Neurons','V2b' = 'Neurons','V3' = 'Neurons','dI6' = 'Neurons')

all <- all[, sample(colnames(all), size = 10000, replace=F)]


##lets downsample it to equal cell type sizes

idents <- Idents(all)
big_clusters <- names(which(table(idents) > 100))  # or some threshold
set.seed(123)
subset_cells <- c()
for (i in levels(idents)) {
  cells_in_cluster <- WhichCells(all, idents = i)
  if (i %in% big_clusters) {
    subset_cells <- c(subset_cells, sample(cells_in_cluster, size = 100))
  } else {
    subset_cells <- c(subset_cells, cells_in_cluster)
  }
}
all <- subset(all, cells = subset_cells)


#prepping resolve spatial data
##We use minmol 6 and scale 30

library(dplyr)
library(data.table)  
library(tidyr)
library(purrr) 


base_dir <- 'DevelopmentSpatialData'

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
  'chat.L' = 'LOC108696296',
  'hoxc6.L' = 'LOC108708639',
  'olig2.L' = 'LOC108707878',
  'otp.L' = 'LOC108717966',
  'pou6f2.L' = 'LOC108719010',
  'sox2.L_10.1' = 'sox2.L'
)

current_rownames <- rownames(spatial)
new_rownames <- current_rownames
new_rownames[current_rownames %in% names(rowname_mapping)] <- rowname_mapping[current_rownames[current_rownames %in% names(rowname_mapping)]]
rownames(spatial) <- new_rownames


##meta data 

  read_and_label <- function(file) {
  dat <- fread(file)

  slide_name <- basename(file)
  slide_name <- sub("_cell_stats\\.csv$", "", slide_name)
  slide_name <- gsub("-", "_", slide_name)

  dat$slide <- slide_name
  dat
}




merged_cell_stats <- map_dfr(cell_stats_files, read_and_label, .id = "source")

merged_cell_stats$source <- NULL


merged_cell_stats$slide <- gsub(".*slide([A-Z0-9\\-_]+)_results_colsadded.*", "\\1", merged_cell_stats$slide)
merged_cell_stats$slide <- gsub("-", "_", merged_cell_stats$slide)

merged_cell_stats$cell <- gsub("-", ".", merged_cell_stats$cell)
merged_cell_stats$slide <- sub(".*_([A-Za-z0-9]+_[0-9]+)$", "\\1", merged_cell_stats$slide)

rownames(merged_cell_stats) <- merged_cell_stats$cell


meta <- merged_cell_stats

#Making Seurat Object with this spatial data

sp <- CreateSeuratObject(spatial, meta.data = meta)


##Now lets remove cells outside gray matter since they are uninformative by drawing borders manually on slides
library(sp)


# Step 1: Initialize an empty vector to store labels for all sections
inout_labels <- rep(NA, nrow(sp))  # Placeholder for the entire dataset

# Step 2: Loop through each section in sp_orig$slide
unique_sections <- unique(sp$slide)  # Get unique section identifiers
for (section in unique_sections) {
  # Subset the data for the current section
  section_indices <- which(sp$slide == section)  # Indices of cells in this section
  section_coords <- data.frame(
    x = sp$x[section_indices],
    y = sp$y[section_indices]
  )
  
  # Step 3: Plot and draw boundary for the current section
  plot(section_coords$y, section_coords$x, main = paste("Draw Boundary for Section:", section))
  boundary_points <- locator(type = "p", pch = 20, col = "red")
  
  # Step 4: Define the boundary polygon
  boundary_coords <- cbind(boundary_points$y, boundary_points$x)  # Flip x and y
  boundary_coords <- rbind(boundary_coords, boundary_coords[1, ])  # Close the boundary
  boundary_polygon <- Polygon(boundary_coords)
  boundary <- Polygons(list(boundary_polygon), "1")
  sp_boundary <- SpatialPolygons(list(boundary))
  
  # Step 5: Determine points inside the boundary
  section_coords_matrix <- cbind(section_coords$x, section_coords$y)  # Original x and y
  inside_boundary <- point.in.polygon(
    section_coords_matrix[, 1],
    section_coords_matrix[, 2],
    boundary_coords[, 1],
    boundary_coords[, 2]
  )
  
  # Step 6: Assign labels ("Inside" or "Outside") for this section
  section_labels <- ifelse(inside_boundary > 0, "Inside", "Outside")
  
  # Step 7: Store these labels in the main inout_labels vector
  inout_labels[section_indices] <- section_labels
}



sp$region <- inout_labels



##Let's order slides so they go by animal and level
Idents(sp) <- sp$slide


sp <- RenameIdents(sp, 'C1_5' = 'C1_5_Animal1_Thoracic', 'C1_6' = 'C1_6_Animal1_Thoracic',
'B2_2' = 'B2_2_Animal2_Brachial', 'B2_3' = 'B2_3_Animal2_Brachial',
'B2_4'= 'B2_4_Animal2_Thoracic','C1_1'= 'C1_1_Animal2_Thoracic','C1_2'= 'C1_2_Animal2_Thoracic',
'C1_3'= 'C1_3_Animal2_Lumbar','C1_4'= 'C1_4_Animal2_Lumbar',
'A1_4'='A1_4_Animal3_Thoracic','A1_5'='A1_5_Animal3_Thoracic','C2_3'='C2_3_Animal3_Thoracic',
'A1_2'= 'A1_2_Animal3_Lumbar','A1_3'= 'A1_3_Animal3_Lumbar','C2_1'= 'C2_1_Animal3_Lumbar')

#lets order by level (brachial, thoracic, lumbar)
levels(sp) <- c('B2_2_Animal2_Brachial', 'B2_3_Animal2_Brachial', 'C1_5_Animal1_Thoracic',
'C1_6_Animal1_Thoracic', 'B2_4_Animal2_Thoracic', 'C1_1_Animal2_Thoracic',
'C1_2_Animal2_Thoracic', 'A1_4_Animal3_Thoracic', 'A1_5_Animal3_Thoracic',
'C2_3_Animal3_Thoracic', 'C1_3_Animal2_Lumbar', 'C1_4_Animal2_Lumbar',
'A1_2_Animal3_Lumbar', 'A1_3_Animal3_Lumbar', 'C2_1_Animal3_Lumbar')

sp$slide <- Idents(sp)


#####Lets rotate slides so they're places the same from dorsal to ventral

manual_angles <- c(
  "B2_2_Animal2_Brachial" = -195, "B2_3_Animal2_Brachial" = -195, "C1_5_Animal1_Thoracic" = 10, "C1_6_Animal1_Thoracic" = 5,
  "B2_4_Animal2_Thoracic" = -195, "C1_1_Animal2_Thoracic" = -195,  "C1_2_Animal2_Thoracic" = -190,  "A1_4_Animal3_Thoracic" = -195,
  "A1_5_Animal3_Thoracic" = 15, "C2_3_Animal3_Thoracic" = 175,  "C1_3_Animal2_Lumbar" = -205,  "C1_4_Animal2_Lumbar" = 170,
  "A1_2_Animal3_Lumbar" = 180,   "A1_3_Animal3_Lumbar" = 180,  "C2_1_Animal3_Lumbar" = 100
)

angles_radians <- manual_angles * pi / 180

# Loop over each slide and apply rotation
for (s in names(angles_radians)) {
  
  # Which cells belong to this slide?
  cells_in_slide <- WhichCells(sp, expression = slide == s)
  
  # Current angle (radians)
  theta <- angles_radians[s]
  
  # Extract x,y
  x <- sp@meta.data[cells_in_slide, "x"]
  y <- sp@meta.data[cells_in_slide, "y"]
  
  # (Optional) center around mean for a more "natural" rotation
  x_center <- mean(x)
  y_center <- mean(y)
  x_adj <- x - x_center
  y_adj <- y - y_center
  
  # 2D rotation transform
  x_rot <- x_adj * cos(theta) - y_adj * sin(theta)
  y_rot <- x_adj * sin(theta) + y_adj * cos(theta)
  
  # Shift back
  x_rot <- x_rot + x_center
  y_rot <- y_rot + y_center
  
  # Store in Seurat metadata (new columns)
  sp@meta.data[cells_in_slide, "x_rot"] <- x_rot
  sp@meta.data[cells_in_slide, "y_rot"] <- y_rot
}


##plotting these
sp$x <- sp$x_rot
sp$y <- sp$y_rot

pdf('Supp1_Spatial_borders_Last.pdf', width = 14, height = 16)

ggplot(sp@meta.data, aes(x = y, y = x, color = region)) +
  geom_point(size = 1, alpha = 0.5) +
  facet_wrap(~ slide) +
  
  scale_color_manual(values = c("Inside" = "tomato", "Outside" = "royalblue")) +

  theme_minimal(base_size = 12) +
  
  theme(
    panel.grid = element_blank(),                      
    panel.background = element_rect(fill = "white"),   
    plot.background = element_rect(fill = "white"),    
    strip.background = element_rect(fill = "white"),   
    strip.text = element_text(color = "black"),        
    legend.position = "right"                          
  )

dev.off()


## and perform label transfer to determine cell_types

DefaultAssay(sp) <- 'RNA'
sp <- NormalizeData(sp, features = rownames(sp))
VariableFeatures(sp) <- rownames(sp)
sp <- ScaleData(sp, features = rownames(sp))
sp <- RunPCA(sp, npcs = 10, features = rownames(sp))

DefaultAssay(all) <- 'RNA'
all <- NormalizeData(all, features = rownames(sp))
VariableFeatures(all) <- rownames(sp)
all <- ScaleData(all, features = rownames(sp))
all <- RunPCA(all, npcs = 10, features = rownames(sp))



anchors <- FindTransferAnchors(reference = all, query = sp, normalization.method = "LogNormalize", features = rownames(sp), dims = 1:10)
predictions.assay <- TransferData(anchorset = anchors, refdata = Idents(all), prediction.assay = TRUE)
sp[['predictions']] <- predictions.assay
sp$prediction.id <- GetTransferPredictions(sp, score.filter = 0.5)

Idents(sp) <- sp$prediction.id
sp <- RenameIdents(sp, 'Differentiating Neurons' = 'Progenitors and Diff', 'Neural Progenitors' = 'Progenitors and Diff')
sp$prediction.id <- Idents(sp)
sp$cell_type <- sp$prediction.id

sp <- subset(sp, idents = c('Unassigned'), invert = TRUE)



sp <- RenameIdents(sp, 'MNs' = 'Neurons')
levels(sp) <- c('RP', 'FP', 'Progenitors and Diff', 'Neurons', 'Oligodendrocytes', 'Other')
sp$prediction.id <- Idents(sp)

##Lets plot results

pdf('Supp2_Spatial_CoarseTypes_Last1.pdf', width = 14, height = 16)
ggplot(sp@meta.data, aes(x = y, y = x, color = prediction.id)) +
  geom_point(size = 1, alpha = 0.65)  +
  facet_wrap(~ slide) +
  scale_color_manual(values = c("Neurons" = "dodgerblue","Progenitors and Diff" = "maroon1", 
  "MNs" = "springgreen4",
  "FP" = "purple", "RP" = "darkblue", "Oligodendrocytes" = "orange", 
  "Other" = "gray87" )) +

  theme_minimal(base_size = 12) +
  
  theme(
    panel.grid = element_blank(),                      
    panel.background = element_rect(fill = "white"),   
    plot.background = element_rect(fill = "white"),    
    strip.background = element_rect(fill = "white"),   
    strip.text = element_text(color = "black"),        
    legend.position = "right"                          
  )

dev.off()



##Now - Neurons - cardinal class label transfer

neu_sp <- subset(sp, idents = 'Neurons')



genes <- c("barhl1.L", "lhx2.S", "lhx9.S", "hmx3.L","isl1.L", "tlx3.L", 
"sall3.L",  "lmx1b.1.S", "wt1.L", "dmrt3.L", "evx1.L", "evx2.S",
 "en1.L", "foxp2.L", "lhx3.L", "vsx2.S", "shox2.S",
 "sox14.L",  "sox21.L", "tal1.S", "gata3.L", "nkx2-2.L", "sim1.S")



neu_sp <- NormalizeData(neu_sp, features = rownames(sp))
VariableFeatures(neu_sp) <- rownames(sp)
neu_sp <- ScaleData(neu_sp, features = rownames(sp))
neu_sp <- RunPCA(neu_sp, npcs = 10, features = rownames(sp))

neu <- subset(neu, idents = 'MNs', invert = TRUE)
neu <- NormalizeData(neu, features = rownames(sp))
VariableFeatures(neu) <- rownames(sp)
neu <- ScaleData(neu, features = rownames(sp))
neu <- RunPCA(neu, npcs = 10, features = rownames(sp))



##lets filter out neurons that don't express at least one cardinal class marker from vector 'genes'
DefaultAssay(neu_sp) <- "RNA" 
counts_mat <- GetAssayData(neu_sp, slot = "counts")
genes.use <- intersect(genes, rownames(counts_mat))
counts_sub <- counts_mat[genes.use, , drop = FALSE]
expr_bin <- counts_sub >= 1  
n_markers_expressed <- Matrix::colSums(expr_bin)
cells_keep <- names(n_markers_expressed)[n_markers_expressed >= 1]
sp_filtered <- subset(neu_sp, cells = cells_keep)


neu_sp <- sp_filtered

DefaultAssay(neu) <- "RNA"
DefaultAssay(neu_sp) <- "RNA"

sc_data_log <- GetAssayData(neu, slot = "data")  
sp_data_log <- GetAssayData(neu_sp, slot = "data") 

genes.use <- genes

genes.use <- intersect(genes.use, rownames(sc_data_log))
genes.use <- intersect(genes.use, rownames(sp_data_log))

sc_data_log <- sc_data_log[genes.use, , drop = FALSE]
sp_data_log <- sp_data_log[genes.use, , drop = FALSE]

sc_clusters <- Idents(neu)
clusters.unique <- levels(sc_clusters)


# for each cluster, compute the mean log-expression across cells in that cluster
cluster_avg_log <- sapply(clusters.unique, function(cl) {
  cells_in_cl <- WhichCells(neu, idents = cl)
  if (length(cells_in_cl) == 0) {
    return(rep(NA, length(genes.use)))
  }
  # mean log expression
  rowMeans(sc_data_log[, cells_in_cl, drop = FALSE], na.rm = TRUE)
})

spot.names <- colnames(sp_data_log)
pred_cluster_corr <- rep(NA, length(spot.names))
names(pred_cluster_corr) <- spot.names

for (i in seq_along(spot.names)) {
  spot_i <- spot.names[i]
  spot_vector <- sp_data_log[, spot_i]
  
  # Skip if all values are NA
  if (all(is.na(spot_vector))) {
    # pred_cluster_corr[spot_i] stays NA
    next
  }
  
  # Compute correlation with each cluster-average vector
  cors <- apply(cluster_avg_log, 2, function(ref_vector) {
    # watch out for zero-variance or all-NA
    if (sd(spot_vector, na.rm = TRUE) == 0 || 
        sd(ref_vector, na.rm = TRUE) == 0) {
      return(NA)
    } else {
      # use pairwise.complete.obs to ignore NA positions
      return(cor(spot_vector, ref_vector, method = "pearson", use = "pairwise.complete.obs"))
    }
  })
  
  # If cors is all NA or empty, we can't pick a best cluster
  if (length(cors) == 0 || all(is.na(cors))) {
    best_cl <- NA
  } else {
    best_cl <- names(cors)[which.max(cors)]
  }
  
  pred_cluster_corr[spot_i] <- best_cl
}

neu_sp$pred_corr <- pred_cluster_corr

saveRDS(neu_sp, 'ResolveDev54_Last_Neurons_CardinalClasses.rds')

mns_sp <- subset(sp, idents = 'MNs')
mns_sp$pred_corr <- 'MNs'

neu_sp <- merge(neu_sp, mns_sp)
Idents(neu_sp) <- neu_sp$pred_corr
neu_sp$neural_types <- neu_sp$pred_corr

##Let's plot results for cardinal classes 
pdf('Supp2_Spatial_CardinalClasses_Last.pdf', width = 14, height = 16)
ggplot(neu_sp@meta.data, aes(x = y, y = x, color = pred_corr)) +
  geom_point(size = 1, alpha = 0.65)  +
  facet_wrap(~ slide) +
  scale_color_manual(values = c("dI1" = "orange",  "dI2" = "gold", "dI3" = "tomato1",  "dI4" = "goldenrod3", "dI5" = "olivedrab2", 
"dI6" = "darkolivegreen", "V0" = "lightslateblue", "V1" = "maroon1",  "V2a" = "violet", "V2b" = "violetred", 
"V3" = "lightpink1", "MNs"  = "green3")) +

  theme_minimal(base_size = 12) +
  
  theme(
    panel.grid = element_blank(),                      
    panel.background = element_rect(fill = "white"),   
    plot.background = element_rect(fill = "white"),    
    strip.background = element_rect(fill = "white"),   
    strip.text = element_text(color = "black"),        
    legend.position = "right"                          
  )

dev.off()


dents(neu_sp) <- neu_sp$pred_corr

##Let's plot evrything else
##Expression of cardinal class markers in neurons

##merge with mns_sp - neu_sp <- merge(neu_sp, mns_sp)

levels(neu_sp) <- rev(c('dI1', 'dI2', 'dI3', 'dI4', 'dI5', 'dI6', 'V0', 'V1', 'V2a', 'V2b', 'V3', 'MNs'))

genes_plot <- c("gad2.L","slc17a7.L", "barhl1.L", "lhx2.S", "hmx3.L","isl1.L", "tlx3.L", 
"sall3.L",  "lmx1b.1.S",  "dmrt3.L", "evx1.L", "evx2.S",
 "en1.L", "foxp2.L", "lhx3.L", "vsx2.S", "shox2.S",
 "sox14.L",  "sox21.L", "tal1.S", "gata3.L", "nkx2-2.L", "sim1.S", 
 "slc18a3.S", "slit2.S") 



DotPlot(neu_sp, features = genes_plot, cols = c('white', 'black'))  + RotatedAxis()
              



