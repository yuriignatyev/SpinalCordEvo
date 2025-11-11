## ===== SIM MATRIX ON INTEGRATED OBJECT  =====

library(ComplexHeatmap)
library(circlize)
library(grid)

knn_graph <- dev@graphs$integrated_nn

## Identify cells by species
Frog_cells  <- which(dev$species == "Frog")
Mouse_cells <- which(dev$species == "Mouse")

## Init similarity matrix
sim_frog_mouse <- matrix(0, nrow = length(neural_types), ncol = length(neural_types))
rownames(sim_frog_mouse) <- paste("Frog",  neural_types)
colnames(sim_frog_mouse) <- paste("Mouse", neural_types)

k.mnn <- 20

# Function to compute the similarity matrix for two species
compute_similarity <- function(cluster1_cells, cluster2_cells, species1, species2) {
  sim_matrix <- matrix(0, nrow = length(neural_types), ncol = length(neural_types))
  rownames(sim_matrix) <- paste(species1, neural_types, sep = " ")
  colnames(sim_matrix) <- paste(species2, neural_types, sep = " ")
  
  for (i in seq_along(neural_types)) {
    cluster1_i <- cluster1_cells[dev$neural_types[cluster1_cells] == neural_types[i]]
    for (j in seq_along(neural_types)) {
      cluster2_j <- cluster2_cells[dev$neural_types[cluster2_cells] == neural_types[j]]
      
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

# Compute the similarity matrix for Frog and Mouse
sim_frog_mouse_dev <- compute_similarity(Frog_cells, Mouse_cells, "Frog", "Mouse")



col_fun <- colorRamp2(
  c(0,   0.02, 0.05, 0.07, 0.10, 0.13, 0.16, 0.20, 0.22),
  c("white","grey97","grey90","grey65","grey55","grey45","grey30","grey15","black")
)




Heatmap(sim_frog_mouse_dev, col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE, 
              show_heatmap_legend = TRUE, row_names_side = "left", 
              name = "Similarity",
              cell_fun = function(j, i, x, y, width, height, fill) {
                if(sim_frog_mouse_dev[i, j] > 0.04)
                  grid.text(sprintf("%.2f", sim_frog_mouse_dev[i, j]), x, y, gp = gpar(fontsize = 8))
              })
