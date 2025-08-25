# Immunogenicity Index & Hot/Cold Classification in Visium
# last updated Aug 25, 2025
# running on R v4.5.0

# --- Setup & helpers
library(Seurat)
library(UCell)
library(FNN)
library(igraph)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

#assumes a processed Visium Seurat object with annotations
#obj <- readRDS("Visium_clean_annot.rds")

DefaultAssay(obj) <- "SCT"

#helper functions
#z-score columns
zcols <- function(df) as.data.frame(scale(df))

#keep only genes present in obj
present <- function(genes, assay = "SCT") {
  allg <- rownames(GetAssayData(obj, assay = assay, slot = "data"))
  intersect(genes, allg)
}

# --- Immunogenic (hot) vs immune-cold signatures
sig_cytolytic <- present(c("GZMA","PRF1"))
sig_ifng      <- present(c("CXCL9","CXCL10","CXCL11","STAT1","IRF1","IDO1"))
sig_mhci      <- present(c("HLA-A","HLA-B","HLA-C","B2M"))
sig_mhcii     <- present(c("HLA-DRA","HLA-DRB1"))
sig_chemokine <- present(c("CCL5","CXCL13","XCL1"))
sig_checkpt   <- present(c("PDCD1","CTLA4","LAG3","TIGIT","ICOS","CD27","CD28"))

sig_tgfb_emt  <- present(c("TGFB1","TGFBR2","COL1A1","COL3A1","VIM","ZEB1"))
sig_wnt       <- present(c("WNT3A","CTNNB1","AXIN2","DKK1"))
sig_hypoxia   <- present(c("VEGFA","HIF1A","CA9","LDHA"))
sig_caf       <- present(c("FAP","PDPN","COL1A1","COL6A3","CXCL12","IL6"))

sets <- list(
  CYT     = sig_cytolytic,
  IFNG    = sig_ifng,
  MHCI    = sig_mhci,
  MHCII   = sig_mhcii,
  CHEMO   = sig_chemokine,
  CHKPT   = sig_checkpt,
  TGFbEMT = sig_tgfb_emt,
  WNT     = sig_wnt,
  HYP     = sig_hypoxia,
  CAF     = sig_caf
)

#score with UCell (robust to library size differences)
obj <- AddModuleScore_UCell(obj, features = sets, name = names(sets))

#aggregate pro- vs anti-immune modules
pro_cols  <- c("CYTCYT", "IFNGIFNG", "MHCIMHCI", "MHCIIMHCII", "CHEMOCHEMO", "CHKPTCHKPT")
anti_cols <- c("TGFbEMTTGFbEMT", "WNTWNT", "HYPHYP", "CAFCAF")

pro  <- zcols(obj@meta.data[, pro_cols, drop = FALSE])
anti <- zcols(obj@meta.data[, anti_cols, drop = FALSE])

obj$ImmunogenicityIndex <- rowMeans(pro, na.rm = TRUE) - rowMeans(anti, na.rm = TRUE)

# --- Spatial smoothing + hot/cold calls
coords <- GetTissueCoordinates(obj)[, c("x","y")]
k <- 8
nn <- get.knn(as.matrix(coords), k = k)$nn.index

smooth_vec <- function(v) {
  out <- v
  for (i in seq_along(v)) out[i] <- mean(c(v[i], v[nn[i,]]), na.rm = TRUE)
  out
}
obj$Immunogenicity_smooth <- smooth_vec(obj$ImmunogenicityIndex)

#quantile thresholds for classification
qs <- quantile(obj$Immunogenicity_smooth, c(.25, .85), na.rm = TRUE)

#hot/intermediate/cold
obj$HotCold <- cut(
  obj$Immunogenicity_smooth,
  breaks = c(-Inf, qs[1], qs[2], Inf),
  labels = c("cold","intermediate","hot")
)
obj$HotCold <- factor(obj$HotCold, levels = c("cold","intermediate","hot"))

#graph-based contiguous patches
el <- cbind(rep(1:nrow(coords), each = k), as.vector(t(nn)))
g  <- simplify(graph_from_edgelist(el, directed = FALSE))

lab_by_mask <- function(mask, prefix) {
  idx <- which(mask)
  sg  <- induced_subgraph(g, vids = idx)
  comp <- components(sg)$membership
  out <- rep(NA_character_, nrow(obj))
  out[idx] <- paste0(prefix, "_", comp)
  out
}
obj$HotPatchID  <- lab_by_mask(obj$HotCold == "hot",  "Hot")
obj$ColdPatchID <- lab_by_mask(obj$HotCold == "cold", "Cold")

# --- Tag malignant states A/B/C
A <- obj@meta.data[,"Malignant cell state A", drop = TRUE]
B <- obj@meta.data[,"Malignant cell state B", drop = TRUE]
C <- obj@meta.data[,"Malignant cell state C", drop = TRUE]
mm <- cbind(A = A, B = B, C = C)
obj$malignant_class <- c("A","B","C")[max.col(mm, ties.method = "first")]

# --- SpaCET immune presence cross-check
DefaultAssay(obj) <- "propMatFromSpaCET"
imm_cols <- intersect(colnames(obj@meta.data), c(
  "T CD8","T CD4","NK","cDC","cDC1 CLEC9A","cDC2 CD1C","cDC3 LAMP3","Macrophage M1"
))
if (length(imm_cols) > 0) {
  obj$ImmunePresence <- rowSums(obj@meta.data[, imm_cols, drop = FALSE], na.rm = TRUE)
}
DefaultAssay(obj) <- "SCT"

# --- Visualization
p1 <- SpatialFeaturePlot(obj, features = "Immunogenicity_smooth", alpha = c(0.1,1)) +
  ggtitle("Immunogenicity Index (smoothed)")

p2 <- SpatialDimPlot(
  obj,
  group.by = "HotCold",
  alpha = c(0.5, 0.5),
  cols = c("cold" = "#1f78b4", "intermediate" = "#bdbdbd", "hot" = "#e31a1c")
) + ggtitle("Hot vs Cold micro-environments")

#class-wise Hot/Cold maps
hotcold_cols <- c(cold = "#2B83BA", intermediate = "#BDBDBD", hot = "#D7191C")
classes <- levels(factor(obj$malignant_class))
plots <- lapply(classes, function(cl) {
  SpatialPlot(
    subset(obj, subset = malignant_class == cl),
    group.by       = "HotCold",
    pt.size.factor = 2,
    alpha          = c(0.95, 1),
    image.alpha    = 0.15,
    stroke         = 0.05,
    cols           = hotcold_cols
  ) +
    ggtitle(paste("Hot/Cold in Malignant", cl)) +
    scale_fill_manual(values = hotcold_cols, limits = names(hotcold_cols), drop = FALSE) +
    scale_color_manual(values = hotcold_cols, limits = names(hotcold_cols), drop = FALSE)
})
p3 <- wrap_plots(plots)

# --- Quantification by malignant class
VlnPlot(obj, features = "Immunogenicity_smooth", group.by = "malignant_class", pt.size = 0) +
  ylab("Immunogenicity (z)") + ggtitle("Index by malignant class")

#fraction table (tidy)
prop_class <- prop.table(table(obj$malignant_class, obj$HotCold), 1)
prop_long <- as.data.frame(as.table(prop_class)) %>%
  rename(malignant_class = Var1, HotCold = Var2, fraction = Freq) %>%
  mutate(
    malignant_class = factor(malignant_class, levels = c("A","B","C")),
    HotCold = factor(HotCold, levels = c("cold","intermediate","hot"))
  )

ggplot(prop_long, aes(x = malignant_class, y = fraction, fill = HotCold)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  scale_fill_manual(values = c("cold" = "#377eb8", "intermediate" = "#cccccc", "hot" = "#e41a1c")) +
  labs(
    x = "Malignant class",
    y = "Fraction of patches",
    fill = "Immune state",
    title = "Distribution of Hot/Cold states across malignant classes"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )


saveRDS(obj,"Visium_clean_scored.rds")