## Summing Xenopus isoforms and mapping to mouse-style genes
# The original script was optimized with tidyverse by Alena Kizenko (former rotation student in Sweeney lab, https://github.com/AlenaKizenko) for better performance 

library('SeuratObject')
library('Seurat')
library('dplyr')
library('ggplot2')
library('tidyverse')
library('readxl')
library('reshape2')
library('RColorBrewer')
library('data.table')

#------------------------------mouse-frog annotation----------------------------

#the gene correspondence table is the current Table S3 of the manuscript
frog_mouse_gene_map = read_xlsx('/.../mouse_frog_gene_correspondence.xlsx',
                                     col_names = FALSE, col_types = 'text') %>%
  select(c(1,2)) %>%
  rename('Xlae_gene' = '...1' , 'Mmus_gene' = '...2') %>%
  as.data.frame()


#---------------------------transforming frog spinal neurons from nf54 -----------
# summing isoforms expression values
frog <- readRDS('/.../Xenopus54_Neurons_CardinalClasses_Final_Labelled.rds')
counts = as.data.frame(GetAssayData(frog, slot = "counts"))
counts = counts[order(rownames(counts)), ]

new_counts = counts %>%
  rownames_to_column('Xlae_gene') %>%
  left_join(frog_mouse_gene_map, by = 'Xlae_gene') %>%
  mutate(Gene = str_to_title(str_remove(string = coalesce(Mmus_gene, Xlae_gene), pattern =  "\\.[L|S]*.*$"))) %>%
  select(-c('Xlae_gene', 'Mmus_gene'))


new_counts = data.table(new_counts)

new_counts = new_counts[ , lapply(.SD, sum), by = Gene]


final_counts = new_counts %>%
  column_to_rownames(var = "Gene")


mfrog = CreateSeuratObject(counts = final_counts)

Idents(mfrog) <- Idents(frog)
mfrog@meta.data <- frog@meta.data

## then the object was downsampled to the mouse neuronal object size (15420 cells) 
mfrog <- mfrog[, sample(colnames(mfrog), size = ncol(15420), replace=F)]


saveRDS(mfrog, file = '/.../Xenopus54_Neurons_CardinalClasses_MouseOrthologs_IntegrationReady.rds')
