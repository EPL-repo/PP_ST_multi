# PP_ST_multi

Companion code for Adimi et al 2024 ("A new spatial multi-omics approach to deeply characterize human cancer tissue using a single tissue section") for the analysis of the generated spatial transcriptomics and proteomics data.

### Data repository
Please refer to the Zenodo repository associated with the manuscript for data access: 10.5281/zenodo.13736222

## Contents

> Directory details

- [Visium ST processing](#STVisium)
- [Spatial alignment of Phenocycler to Visium data](#ImageAlignment)
- [Importing processed Phenocycler data into R]

---

## Visium ST processing

---

## Spatial alignment of Phenocycler to Visium data

---

## Importing processed Phenocycler data into R

- This can be performed using the [Giotto R package](https://giottosuite.readthedocs.io/en/master/)

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
