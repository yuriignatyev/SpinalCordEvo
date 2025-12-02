Code for 'Innovations in Spinal Cord Cell Type Heterogeneity Across Vertebrate Evolution' by *Ignatyev et al.*, 2025 (https://www.biorxiv.org/content/10.1101/2025.10.09.680955v2)

Datasets are currently available via private link and will be made available upon publication. The repository will be updated with additional scripts as needed. 

---

## 📘 Overview
This repository contains the analysis scripts and notebooks used to produce the results and figures in the manuscript.  
The study integrates developmental and adult single-cell/nucleus and spatial transcriptomic data from **Xenopus laevis**, **mouse**, **human**, and **zebrafish** to explore the evolutionary origins of spinal cord neuronal diversity.

![](abstract.jpg)

---

## 📁 Structure
The repository is organized into subdirectories corresponding to the analyses used to generate the figures:

- **scripts/**
  - Contains codes used for data processing and cross-species analyses.
  - Analyses are grouped into:
    - **Development** (Figures **1–2**, Supplementary **S1–S9**) <br>
      • data processing (01.-03.) <br>
      • transforming frog data into mouse format (04. - was also used for transforming **adult** frog data) <br>
      • purity analysis and Pearson correlation of cardinal class expression profiles (05.) - Fig.2D, Supp4 <br>
      • developmental integration (06.) - Fig.2E, Supp3 <br>
      • similarity index in development (07.), Fig.2E, Supp3 <br>
      • frog spatial data analysis (08.), Fig.2F-I, Supp7-9 <br>
      
    - **Adult** (Figures **3–4**, Supplementary **S11–S18**) <br>
      • frog adult snrna-seq data processing  (01.) <br> -  Fig.3
      • frog adult spatial data processing (02.) <br> - Fig.3
      • adult integration (03.) <br>
      • similarity index in adult (04.) - Fig.4F <br>
      • species mixing (05.), Fig.4, Supp18 <br>
      • projecting divergent and conserved populations onto frog and mouse spatial data (06.), Fig.4J
    
---

## Dependencies

The following core R packages were used:

- **Seurat**
- **dplyr**
- **ggplot2**
- **tidyr**
- **data.table**
- **Matrix**
- **patchwork**
- **cowplot**
- **mixtools**
- **tibble**
- **biomaRt**
- **stringr**
- **grid**
- **purrr**
- **sp**
- **ggtern**
- **ggrepel**
- **viridis**
- **RColorBrewer**
- **gridExtra**
- **future**
- **future.apply**
- **ComplexHeatmap**
- **CellChat**
- **circlize**



