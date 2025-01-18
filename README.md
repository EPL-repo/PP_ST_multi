# PP_ST_multi

Companion code of Adimi et al 2024 **"A new spatial multi-omics approach to deeply characterize human cancer tissue using a single tissue section"**, for the analysis of the generated spatial transcriptomics and proteomics data. Preprint can be found [here](https://www.biorxiv.org/content/10.1101/2024.10.03.616510v1).

### Data repository
Please refer to the Zenodo repository associated with the manuscript for data access: [10.5281/zenodo.13736222](https://zenodo.org/records/13736222)

## Contents

> Directory details

- [Visium ST processing](#STVisium)
- [Spatial alignment of Phenocycler to Visium data](#ImageAlignment)
- Importing processed Phenocycler data into R

---

## Visium ST processing
- Two Jupyter notebooks provided to detail Visium SD data processing (*STVisium_preprocessing_manuscript*) and spot deconvolution using SpaCET (*STVisium_Malignant_deconvolution_manuscript*)
- Both routines running smoothly on:
```shell
R version 4.4.1 (2024-06-14)
Platform: aarch64-apple-darwin20
Running under: macOS Sonoma 14.6.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Europe/Paris
tzcode source: internal 
```
---

## Spatial alignment of Phenocycler to Visium data
- The image alignment and generation of a pseudo-spot grid for Phenocycler data procedures were performed in **Python** (v3.11.9).
- The alignment procedure involved identifying the optimal affine transformation using two parameters: 1) the scaling factor and 2) the shift term. The scaling factor was determined either directly from the microscopy images’ resolution metadata or estimated from the images themselves. The estimation protocol involved comparing the distribution of vertical and horizontal distances between the masked Visium and Phenocycler images. The shift term was determined by aligning the centers of mass of the two masked images.
These parameters were then fine-tuned using a grid search optimization procedure. The masked images were aligned, and the matching score was calculated using one of two alignment scores: 1) the total number of mismatched pixels in the aligned masks, or 2) the geometric mean of the number of mismatched pixels from each of the two masks independently.
- For the creation of a pseudo-spot grid for Phenocycler data, the identified affine image transformation was applied to the cell centroids segmented from the Phenocycler data. The transformed centroid coordinates were then assigned to the coordinates of Visium spots by matching each centroid to the nearest spot within the pseudospot boundary (defined as a distance of less than 27.5 μm from the cell centroid to the spot center). The cell-level Phenocycler expression signals were then aggregated to the pseudo-spot level by summing the signals from all the corresponding cells. The identity of a pseudospot was defined by the most abundant cell type or state within each pseudospot.

---

## Importing processed Phenocycler data into R
- Data generated through the Enable Medicine Portal / ATLAS can be accessed and exported using Enable's Workbench and their two R packages: emconnect and SpatialMap.
- Using SpatialMap, we exported structured data to CSV files containing: segmented cell centroids ("coords.csv"), expression data ("expr.csv") and all the generated metadata - such as annotations, clusters, etc. ("meta.csv") for each high-plex image
- Importing these for further processing can be performed using the [Giotto R package](https://giottosuite.readthedocs.io/en/master/)

```shell
#--- Setup
library(Giotto)

results_folder = "~/PP_ST_multi/codex_eg/"

#python_path = NULL 
#if(is.null(python_path)) {
#  installGiottoEnvironment()
#}

#--- Global instructions
instrs = createGiottoInstructions(show_plot = FALSE,
                                  save_plot = TRUE,
                                  save_dir = results_folder)
                                
expr_path = paste0(results_folder, "PP/expr.csv")
loc_path  = paste0(results_folder, "PP/coords.csv")
meta_path = paste0(results_folder, "PP/meta.csv")

#tile size to adjust coordinates
xtilespan = 1344;
ytilespan = 1008;

#--- Create Giotto object and pre-process
codex_expression = readExprMatrix(expr_path, transpose = T)
codex_locations  = data.table::fread(loc_path)
codex_metadata   = data.table::fread(meta_path)

stitch_file = stitchTileCoordinates(location_file = codex_metadata, Xtilespan = xtilespan, Ytilespan = ytilespan);
codex_locations = stitch_file[,.(Xcoord, Ycoord)]

codex_test <- createGiottoObject(raw_exprs = codex_expression, 
                                 spatial_locs = codex_locations,
                                 instructions = instrs,
                                 cell_metadata = codex_metadata)

#An object of class giotto 
#>Active spat_unit:  cell 
#>Active feat_type:  rna 
#[SUBCELLULAR INFO]
#[AGGREGATE INFO]
#expression -----------------------
#  [cell][rna] raw
#spatial locations ----------------
#  [cell] raw

## visualize
spatPlot(gobject = codex_test,point_size = 0.1, 
         coord_fix_ratio = NULL,point_shape = 'no_border',
         save_param = list(save_name = 'spatPlot_DATA'))
                  
## filter
DG_subset <- filterGiotto(gobject = codex_test,
                           expression_threshold = 1,
                           gene_det_in_min_cells = 10,
                           min_det_genes_per_cell = 2,
                           expression_values = c('raw'),
                           verbose = T)

DG_subset <- normalizeGiotto(gobject = codex_test, scalefactor = 6000, verbose = T,
                              log_norm = FALSE,library_size_norm = FALSE,
                              scale_feats = FALSE, scale_cells = TRUE)

## add gene & cell statistics
codex_test <- addStatistics(gobject = codex_test,expression_values = "normalized")


```


---
