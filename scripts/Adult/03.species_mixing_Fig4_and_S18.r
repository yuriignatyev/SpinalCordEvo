library(Seurat)
library(dplyr)
library(ggplot2)
library(ggtern)
library(patchwork) 

# Int here is your cross-species integrated object, for example, dorsal excitatory integrated neurons subset
int <- readRDS('/.../DorsalExcitatoryNeurons_Integrated_Subset_Across_FrogMouseHuman_CCA.rds')

int$species <- as.factor(int$species)

# 1. Calculate Species Counts per Cluster
cluster_species_counts <- table(int$seurat_clusters, int$species)

# 2. Convert to Data Frame
cluster_species_df <- as.data.frame.matrix(cluster_species_counts) %>%
  mutate(cluster = rownames(.))

# 3. Normalize Counts by Total Species Counts (Proportional Normalization)
total_species <- colSums(cluster_species_counts)

cluster_species_df <- cluster_species_df %>%
  mutate(
    Frog_prop = Frog / total_species["Frog"],
    Human_prop = Human / total_species["Human"],
    Mouse_prop = Mouse / total_species["Mouse"]
  )

# 4. Calculate Species Proportions Within Each Cluster
cluster_species_df <- cluster_species_df %>%
  rowwise() %>%
  mutate(
    total_norm = Frog_prop + Human_prop + Mouse_prop,
    frog_prop = ifelse(total_norm > 0, Frog_prop / total_norm, 0),
    human_prop = ifelse(total_norm > 0, Human_prop / total_norm, 0),
    mouse_prop = ifelse(total_norm > 0, Mouse_prop / total_norm, 0)
  ) %>%
  ungroup()



# 5. Apply Square Root Transformation to Species Proportions
cluster_species_df <- cluster_species_df %>%
  mutate(
    frog_prop_sqrt = sqrt(frog_prop),
    human_prop_sqrt = sqrt(human_prop),
    mouse_prop_sqrt = sqrt(mouse_prop)
  )

# 6. Normalize the Transformed Proportions
cluster_species_df <- cluster_species_df %>%
  rowwise() %>%
  mutate(
    total_prop_sqrt = frog_prop_sqrt + human_prop_sqrt + mouse_prop_sqrt,
    frog_prop_adj = ifelse(total_prop_sqrt > 0, frog_prop_sqrt / total_prop_sqrt, 0),
    human_prop_adj = ifelse(total_prop_sqrt > 0, human_prop_sqrt / total_prop_sqrt, 0),
    mouse_prop_adj = ifelse(total_prop_sqrt > 0, mouse_prop_sqrt / total_prop_sqrt, 0)
  ) %>%
  ungroup()

# 7. Define Base Colors for Each Species (RGB, [0,1] scale)
frog_color  <- c(0, 255, 0) / 255      # Green
human_color <- c(255, 0, 0) / 255      # Red
mouse_color <- c(0, 0, 255) / 255      # Blue

# 8. Combine Adjusted Proportions with Base Colors
cluster_species_df <- cluster_species_df %>%
  rowwise() %>%
  mutate(
    R = human_prop_adj * human_color[1],  # Human contributes to Red
    G = frog_prop_adj * frog_color[2],    # Frog contributes to Green
    B = mouse_prop_adj * mouse_color[3]    # Mouse contributes to Blue
  ) %>%
  ungroup()

## 9. Calculate Standard Deviation of Adjusted Proportions
cluster_species_df <- cluster_species_df %>%
  rowwise() %>%
  mutate(
    sd_prop = sd(c(frog_prop_adj, human_prop_adj, mouse_prop_adj))
  ) %>%
  ungroup()

# 10. Define Blend Factor Based on Minimum Species Proportion
# Define thresholds
threshold_min <- 0.05  
threshold_max <- 0.20  

blend_power <- 0.5  # Controls the smoothness of blending

cluster_species_df <- cluster_species_df %>%
  mutate(
    # Calculate the minimum species proportion
    min_prop = pmin(frog_prop, human_prop, mouse_prop),
    
    # Calculate blend_factor based on min_prop and thresholds
    blend_factor = case_when(
      min_prop <= threshold_min ~ 0,
      min_prop >= threshold_max ~ 1,
      TRUE ~ ((min_prop - threshold_min) / (threshold_max - threshold_min)) ^ blend_power
    ),
    
    # Ensure blend_factor is within [0,1]
    blend_factor = pmin(pmax(blend_factor, 0), 1)
  )

# 11. Blend RGB Values Towards Gray
gray_point <- 0.75  # Adjusted gray shade

cluster_species_df <- cluster_species_df %>%
  mutate(
    R_adj = R * (1 - blend_factor) + blend_factor * gray_point,
    G_adj = G * (1 - blend_factor) + blend_factor * gray_point,
    B_adj = B * (1 - blend_factor) + blend_factor * gray_point,
    # Ensure RGB values are within [0,1], set to gray_point if NaN
    R_adj = ifelse(is.nan(R_adj) | is.na(R_adj), gray_point, pmin(pmax(R_adj, 0), 1)),
    G_adj = ifelse(is.nan(G_adj) | is.na(G_adj), gray_point, pmin(pmax(G_adj, 0), 1)),
    B_adj = ifelse(is.nan(B_adj) | is.na(B_adj), gray_point, pmin(pmax(B_adj, 0), 1)),
    color = rgb(R_adj, G_adj, B_adj)
  )

# 12. Assign Colors to Clusters
cluster_colors <- cluster_species_df$color
names(cluster_colors) <- cluster_species_df$cluster

# 13. Create a Vector of Colors for Each Cell Based on Cluster
cell_cluster_colors <- cluster_colors[as.character(int$seurat_clusters)]
names(cell_cluster_colors) <- colnames(int)

# 14. Add Cluster Colors to Seurat Object's Metadata
int$cluster_color <- cell_cluster_colors

# 15. Extract UMAP Embeddings
umap_df <- Embeddings(int, "umap") %>% as.data.frame()

# Ensure the row names match between embeddings and metadata
if (!identical(rownames(umap_df), colnames(int))) {
  umap_df <- umap_df[match(colnames(int), rownames(umap_df)), ]
}

# 16. Add Cluster Colors and Cluster IDs to UMAP Data Frame
umap_df <- umap_df %>%
  mutate(
    cluster_color = int$cluster_color,
    seurat_clusters = int$seurat_clusters
  )

# Rename UMAP Columns if Necessary
colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")

# 17. Generate UMAP Plot with Custom Color Scheme
umap_plot <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = cluster_color)) +
  geom_point(size = 0.5) +
  scale_color_identity() +
  theme_void() +
  labs(title = "UMAP of Clusters Colored by Species Mixing")

# 18. Create Ternary Plot as Legend
# Prepare data for ternary plot
ternary_df <- cluster_species_df %>%
  select(cluster, frog_prop, human_prop, mouse_prop, color)

# Generate a grid for the ternary gradient
grid_size <- 200  # Increased grid size for smoother gradient
grid_df <- expand.grid(
  frog_prop = seq(0, 1, length.out = grid_size),
  human_prop = seq(0, 1, length.out = grid_size)
) %>%
  # Ensure that frog_prop + human_prop + mouse_prop = 1
  mutate(mouse_prop = 1 - frog_prop - human_prop) %>%
  filter(mouse_prop >= 0)  # Keep valid proportions

# Apply the same color assignment logic to the grid
grid_df <- grid_df %>%
  rowwise() %>%
  mutate(
    # Calculate the minimum species proportion
    min_prop = pmin(frog_prop, human_prop, mouse_prop),
    
    # Calculate blend_factor based on min_prop and thresholds
    blend_factor = case_when(
      min_prop <= threshold_min ~ 0,
      min_prop >= threshold_max ~ 1,
      TRUE ~ ((min_prop - threshold_min) / (threshold_max - threshold_min)) ^ blend_power
    ),
    
    # Ensure blend_factor is within [0,1]
    blend_factor = pmin(pmax(blend_factor, 0), 1),
    
    # Calculate RGB based on species proportions
    R = human_prop * human_color[1],  # Human contributes to Red
    G = frog_prop * frog_color[2],    # Frog contributes to Green
    B = mouse_prop * mouse_color[3],  # Mouse contributes to Blue
    
    # Blend towards gray
    R_adj = R * (1 - blend_factor) + blend_factor * gray_point,
    G_adj = G * (1 - blend_factor) + blend_factor * gray_point,
    B_adj = B * (1 - blend_factor) + blend_factor * gray_point,
    
    # Ensure RGB values are within [0,1], set to gray_point if NaN
    R_adj = ifelse(is.nan(R_adj) | is.na(R_adj), gray_point, pmin(pmax(R_adj, 0), 1)),
    G_adj = ifelse(is.nan(G_adj) | is.na(G_adj), gray_point, pmin(pmax(G_adj, 0), 1)),
    B_adj = ifelse(is.nan(B_adj) | is.na(B_adj), gray_point, pmin(pmax(B_adj, 0), 1)),
    
    # Assign color based on blend_factor
    color = rgb(R_adj, G_adj, B_adj)
  ) %>%
  ungroup()

# Generate the ternary gradient plot
ternary_gradient <- ggtern(data = grid_df, aes(x = frog_prop, y = human_prop, z = mouse_prop)) +
  geom_point(aes(color = color), size = 1, alpha = 0.3) +  # Background gradient
  scale_color_identity() +
  theme_bw() +
  theme_showarrows() +
  labs(title = "Triangular Gradient Legend",
       T = "Frog",
       L = "Human",
       R = "Mouse") +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5)
  ) +
  # Overlay cluster points as black dots
  geom_point(data = ternary_df, 
             aes(x = frog_prop, y = human_prop, z = mouse_prop), 
             color = "black", 
             size = 2) +
  # Add cluster number labels on top of dots
  geom_text(data = ternary_df, 
            aes(x = frog_prop, y = human_prop, z = mouse_prop, label = cluster), 
            color = "black", 
            size = 3, 
            vjust = -1, 
            hjust = 0.5)

# 20. Combine UMAP Plot and Ternary Gradient Legend
combined_plot <- umap_plot + ternary_gradient +
  plot_layout(ncol = 2, widths = c(3, 2)) +  # Adjust widths as needed
  plot_annotation(title = "UMAP with Triangular Gradient Legend",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 16)))

# 21. Display the Combined Plot
combined_plot 
