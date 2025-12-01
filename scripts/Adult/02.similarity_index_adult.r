library(Seurat)

#loading interated neuronal object
int <- readRDS('/.../Neurons_Integrated_Across_FrogMouseHuman_CCA.rds')

int$neural_types <- int$coarse
int <- FindNeighbors(int, dims = 1:150, k.param = 20)

neural_types <- c('Exc', 'Inh', 'Vent')
k.mnn <- 20

knn_graph <- int@graphs$integrated_nn

# Identify cells by species
Frog_cells <- which(int$species == "Frog")
Mouse_cells <- which(int$species == "Mouse")
Human_cells <- which(int$species == "Human")

# Initialize the similarity matrices for each pairwise comparison
sim_frog_mouse <- matrix(0, nrow = length(neural_types), ncol = length(neural_types))
sim_frog_human <- matrix(0, nrow = length(neural_types), ncol = length(neural_types))
sim_human_mouse <- matrix(0, nrow = length(neural_types), ncol = length(neural_types))

rownames(sim_frog_mouse) <- colnames(sim_frog_mouse) <- paste("Frog", neural_types, sep = " ")
rownames(sim_frog_human) <- colnames(sim_frog_human) <- paste("Frog", neural_types, sep = " ")
rownames(sim_human_mouse) <- colnames(sim_human_mouse) <- paste("Human", neural_types, sep = " ")

colnames(sim_frog_mouse) <- paste("Mouse", neural_types, sep = " ")
colnames(sim_frog_human) <- paste("Human", neural_types, sep = " ")
rownames(sim_human_mouse) <- paste("Mouse", neural_types, sep = " ")

compute_similarity <- function(cluster1_cells, cluster2_cells, species1, species2) {
  sim_matrix <- matrix(0, nrow = length(neural_types), ncol = length(neural_types))
  rownames(sim_matrix) <- paste(species1, neural_types, sep = " ")
  colnames(sim_matrix) <- paste(species2, neural_types, sep = " ")
  
  for (i in seq_along(neural_types)) {
    cluster1_i <- cluster1_cells[int$neural_types[cluster1_cells] == neural_types[i]]
    for (j in seq_along(neural_types)) {
      cluster2_j <- cluster2_cells[int$neural_types[cluster2_cells] == neural_types[j]]
      
      if (length(cluster1_i) > 0 && length(cluster2_j) > 0) {
        # Compute the overlap for cluster1 to cluster2 neural types
        cluster1_to_cluster2_overlap <- rowSums(knn_graph[cluster1_i, cluster2_j] > 0) / k.mnn
        cluster1_to_cluster2_avg <- mean(cluster1_to_cluster2_overlap)
        
        # Compute the overlap for cluster2 to cluster1 neural types
        cluster2_to_cluster1_overlap <- rowSums(knn_graph[cluster2_j, cluster1_i] > 0) / k.mnn
        cluster2_to_cluster1_avg <- mean(cluster2_to_cluster1_overlap)
        
        # Average the overlaps to get a similarity score for the (i, j) pair
        sim_matrix[i, j] <- (cluster1_to_cluster2_avg + cluster2_to_cluster1_avg) / 2
      }
    }
  }
  return(sim_matrix)
}

# Compute similarity matrices for each pairwise comparison

sim_frog_mouse <- compute_similarity(Frog_cells, Mouse_cells, "Frog", "Mouse")
sim_frog_human <- compute_similarity(Frog_cells, Human_cells, "Frog", "Human")
sim_human_mouse <- compute_similarity(Human_cells, Mouse_cells, "Human", "Mouse")



##Heatmap
library(ComplexHeatmap)
library(CellChat)
library(circlize)

color.heatmap <- colorRamp3(c(0.01 , 0.03, 0.07, 0.1, 0.13, 0.15, 0.23), c('grey100', 'grey95','grey90','grey80','grey26','grey10','grey0'))


Heatmap(sim_frog_mouse, col = color.heatmap, cluster_rows = FALSE, cluster_columns = FALSE, 
              show_heatmap_legend = TRUE, row_names_side = "left", 
              name = "Similarity",
              cell_fun = function(j, i, x, y, width, height, fill) {
                if(sim_frog_mouse[i, j] > 0.04)
                  grid.text(sprintf("%.2f", sim_frog_mouse[i, j]), x, y, gp = gpar(fontsize = 8))
              })



Heatmap(sim_frog_human, col = color.heatmap, cluster_rows = FALSE, cluster_columns = FALSE, 
              show_heatmap_legend = TRUE, row_names_side = "left", 
              name = "Similarity",
              cell_fun = function(j, i, x, y, width, height, fill) {
                if(sim_frog_human[i, j] > 0.04)
                  grid.text(sprintf("%.2f", sim_frog_human[i, j]), x, y, gp = gpar(fontsize = 8))
              })

Heatmap(sim_human_mouse, col = color.heatmap, cluster_rows = FALSE, cluster_columns = FALSE, 
              show_heatmap_legend = TRUE, row_names_side = "left", 
              name = "Similarity",
              cell_fun = function(j, i, x, y, width, height, fill) {
                if(sim_human_mouse[i, j] > 0.00)
                  grid.text(sprintf("%.2f", sim_human_mouse[i, j]), x, y, gp = gpar(fontsize = 8))
              })


