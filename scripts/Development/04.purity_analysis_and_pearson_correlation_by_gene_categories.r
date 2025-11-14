library(biomaRt)
library(Seurat)
library(stringr)
library(ggplot2)

#for individual object, whether mouse or mouse-transformed frog
ensembl1 <- useMart("ensembl",dataset="hsapiens_gene_ensembl", host = "useast.ensembl.org")
x <- getBM(attributes=c('hgnc_symbol', 'go_id', 'name_1006', "namespace_1003"), filters = 'hgnc_symbol', values = genes, mart= ensembl1)

##TFs
y <- x[grep("DNA-binding transcription", x$name_1006),]
TFs_list1 <- unique(y$hgnc_symbol)
y <- x[grep("regulation of transcription, DNA-templated", x$name_1006),]
TFs_list2 <- unique(y$hgnc_symbol)
list12 <- c(TFs_list1, TFs_list2)
unique <- unique(list12)
TFs_list <- str_to_title(unique)
##Ion Transport
y <- x[grep(c("ion transport"), x$name_1006),]
ion_list <- str_to_title(unique(y$hgnc_symbol))
#singaling
y <- x[grep(c("signaling receptor"), x$name_1006),]
signaling_receptor_list <- str_to_title(unique(y$hgnc_symbol))
#CA
y <- x[grep(c("cell adhesion"), x$name_1006),]
cell_adhesion_list <- str_to_title(unique(y$hgnc_symbol))



####PURITY

categories <- list(
  TFs = TFs_list,
  cell_adhesion = cell_adhesion_list,
  ion_transport = ion_list,
  signaling = signaling_receptor_list,
  allgenes = VariableFeatures(seu)
)

# Prepare a data frame for purity results
purities <- data.frame(cluster = character(), purity = numeric(), category = character())

# Process each category
for (cat_name in names(categories)) {
  genes <- categories[[cat_name]]
  
  # Subset Seurat object by genes in the current category
  sub_data <- RunPCA(seu, features = genes, verbose = FALSE)
 
  
  # Perform clustering - adjust parameters as necessary
  sub_data <- FindNeighbors(sub_data, verbose = FALSE)
  sub_data <- FindClusters(sub_data, resolution = 7, verbose = FALSE)
  
  # Calculate purity
  for (i in levels(sub_data$seurat_clusters)) {
    cluster_cells <- WhichCells(sub_data, expression = seurat_clusters == i)
    cell_types <- sub_data$neural_types[cluster_cells]  # assuming 'cell_type' is a defined metadata
    freq_table <- table(cell_types)
    max_freq <- max(freq_table)
    total_cells <- length(cluster_cells)
    purity <- max_freq / total_cells
    
    # Append results
    purities <- rbind(purities, data.frame(cluster = i, purity = purity, category = cat_name))
  }
}


#plotting 
ggplot(purities, aes(x = purity, fill = category)) +
  geom_density(alpha = 0.3, adjust = 1.5) +  # Adjust transparency with alpha and smoothing with adjust
  scale_fill_manual(values = c("deeppink", "green", "blue", "yellow", "tomato")) +  # Customize colors
  xlim(0, 1) +
  theme_minimal() +
  theme(
    text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(hjust = 1, size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )


###PEARSON CORRELATIONS ACROSS SPECIES

##first let's make sure to overlap (intersect) gene lists
##also, since we use scale.data slot, lets scale both of the data objects by their rownames in 'counts' slot 

DefaultAssay(mouse) <- 'RNA'
length(rownames(mouse))

genes <- intersect(rownames(mfrog), rownames(mouse))

DefaultAssay(mouse) <- 'integrated'

mfrog <- ScaleData(mfrog, features = genes)
mouse <- ScaleData(mouse, features = genes)



classes <- c('dI1', 'dI2', 'dI3', 'dI4', 'dI5', 'dI6', 'V0', 'V1', 'V2a', 'V2b', 'V3', 'MNs')

##frog; DefaultAssay(frog) <- 'RNA'
p <- list()
m <- list()
l <- list()

for (i in classes){
  
  p[[i]] <- subset(mfrog, idents = c(i))
  m[[i]] <- as.data.frame(GetAssayData(p[[i]], slot = "scale.data"), genes = rownames(GetAssayData(mfrog, slot = "scale.data")), fix_names = TRUE)
  l[[i]] <- as.data.frame(rowMeans(m[[i]][TFs_list,]))
}


x <- list()
y <- list()
z <- list()

for (i in classes){

  x[[i]] <- subset(mouse, idents = c(i))
  y[[i]] <- as.data.frame(GetAssayData(x[[i]], slot = "scale.data"), genes = rownames(GetAssayData(mouse, slot = "scale.data")), fix_names = TRUE)
  z[[i]] <- as.data.frame(rowMeans(y[[i]][TFs_list,]))
  
  }


M <- cbind(l$dI1, l$dI2, l$dI3, l$dI4, l$dI5, l$dI6, l$V0, l$V1, l$V2a, l$V2b, l$V3, l$MNs,
           z$dI1, z$dI2, z$dI3, z$dI4, z$dI5, z$dI6, z$V0, z$V1, z$V2a, z$V2b, z$V3, z$MNs)

frogclasses <- paste0(classes, "_F")
mouseclasses <- paste0(classes, "_M")

colnames(M) <- append(frogclasses, mouseclasses) 

Pearson <- cor(M, method = "pearson")
Pearson <- Pearson[frogclasses, mouseclasses]


###Plotting

color.heatmap <- colorRamp3(c(0, 0.1,0.2, 0.30,0.45,0.6), c('white', 'gray90','gray80','gray70','gray60','black'))

Heatmap(Pearson,col = color.heatmap, row_names_gp = gpar(fontsize = 18), 
        column_names_gp = gpar(fontsize = 18),cluster_rows = FALSE, cluster_columns = FALSE,show_heatmap_legend = T, row_names_side = "left",column_title = "Pearson", name = NULL,
             cell_fun = function(j, i, x, y, width, height, fill) {
               if(Pearson[i, j] > 0.04)
                 grid.text(sprintf("%.2f", Pearson[i, j]), x, y, gp = gpar(fontsize = 15, col = 'white'))
             })



##default_colors
Heatmap(Pearson, row_names_gp = gpar(fontsize = 18), 
        column_names_gp = gpar(fontsize = 18),cluster_rows = FALSE, cluster_columns = FALSE,show_heatmap_legend = T, row_names_side = "left",column_title = "Pearson", name = NULL,
             cell_fun = function(j, i, x, y, width, height, fill) {
               if(Pearson[i, j] > 0.04)
                 grid.text(sprintf("%.2f", Pearson[i, j]), x, y, gp = gpar(fontsize = 15, col = 'white'))
             })


color.heatmap <- colorRamp3(c(-0.4097346, 0.6213057), c('white', 'black'))

