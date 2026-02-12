library(FNN)
library(dplyr)
library(Seurat)
jsd <- function(p, q, base = exp(1)) {
  # p, q: numeric vectors, non-negative, sum to 1
  p <- p / sum(p)
  q <- q / sum(q)
  m <- 0.5 * (p + q)
  
  kl <- function(a, b) sum(ifelse(a > 0, a * log(a / b, base = base), 0))
  0.5 * kl(p, m) + 0.5 * kl(q, m)
}

profile_to_prob <- function(x, eps = 1e-6) {
  x <- pmax(x, 0)
  x <- x + eps
  x / sum(x)
}

build_knn_rings <- function(
    dd,
    nrings = 10,
    source_label = "TrabBone",
    allowed_label = "Cavity",
    ring_prefix = "Ring"
){
  # Targets (source) and allowed-to-assign spots
  target_spots  <- rownames(dd@meta.data)[dd@meta.data$spa == source_label]
  allowed_spots <- rownames(dd@meta.data)[dd@meta.data$spa == allowed_label]
  
  coords <- GetTissueCoordinates(dd)
  coords_mat <- as.matrix(coords[, c("x", "y")])
  spot_names <- rownames(coords_mat)
  
  target_idx <- match(target_spots, spot_names)
  target_idx <- target_idx[!is.na(target_idx)]
  if (length(target_idx) == 0) stop("No target spots found for source_label = ", source_label)
  
  k_values <- k_hex_cumulative(nrings)
  max_k <- max(k_values$k)
  
  # Compute KNN once
  nn <- get.knn(coords_mat, k = max_k)
  
  assigned_idx <- integer(0)  # track spots already assigned to a ring
  
  for (i in seq_len(nrow(k_values))) {
    r  <- k_values$ring[i]
    kr <- k_values$k[i]
    kprev <- if (i == 1) 0 else k_values$k[i - 1]
    
    # cumulative neighbor indices up to kr / kprev
    neigh_r    <- unique(as.vector(nn$nn.index[target_idx, seq_len(kr), drop = FALSE]))
    neigh_prev <- if (kprev == 0) integer(0) else unique(as.vector(nn$nn.index[target_idx, seq_len(kprev), drop = FALSE]))
    
    # true ring = new shell only
    ring_idx <- setdiff(neigh_r, neigh_prev)
    
    # remove targets themselves (optional safety) and keep only allowed tissue
    ring_spots <- setdiff(spot_names[ring_idx], target_spots)
    ring_spots <- intersect(ring_spots, allowed_spots)
    
    # assign only those not already assigned (closest ring wins)
    ring_spots <- setdiff(ring_spots, spot_names[assigned_idx])
    if (length(ring_spots) == 0) next
    
    dd$spa[ring_spots] <- paste0(ring_prefix, r)
    assigned_idx <- c(assigned_idx, match(ring_spots, spot_names))
  }
  
  dd$idents <- paste(dd$spa, dd$group)
  dd
}
##### Load Data #####

dd <- Load10X_Spatial(data.dir = path.in)
dd@images$slice1@scale.factors$spot <- 6/(dd@images$slice1@scale.factors$lowres) #correct image scale
spa <- read.csv(path.spa, stringsAsFactors = F, header = T, row.names = 1) #load metadata
dd <- AddMetaData(dd, spa) # add metadata
dd <- dd %>% subset(!(spa %in% c('out', 'Undefined', 'CortBone', 'GP', 'TrabBone2')) &
                      nCount_Spatial > 250)
SpatialFeaturePlot(dd, features = 'nFeature_Spatial', max.cutoff = 3000)

#Sample-wise SCT normalization
dd.norm <- SCT2_norm(dd)


dd.norm <- build_knn_rings(dd.norm, nrings = 15)
SpatialDimPlot(dd.norm, group.by = 'spa')
seu_split <- SplitObject(dd.norm, split.by = 'group')
seu_wt <- seu_split$Wt
seu_oe <- seu_split$WntOverExp


library(Seurat)
library(Matrix)

ring_jsd_seurat <- function(seu_wt, seu_oe,
                            assay = DefaultAssay(seu_wt),
                            layer = "data",          # log-normalized
                            ring_col = "spa",
                            min_detect_frac = 0.05, # per section
                            eps = 1e-6) {
  
  X_wt <- GetAssayData(seu_wt, assay = assay, layer = layer)
  X_oe <- GetAssayData(seu_oe, assay = assay, layer = layer)
  
  rings_wt <- seu_wt@meta.data[[ring_col]]
  rings_oe <- seu_oe@meta.data[[ring_col]]
  
  rings_wt <- as.factor(rings_wt)
  rings_oe <- as.factor(rings_oe)
  
  # Ensure same ring levels in both
  all_levels <- sort(unique(c(levels(rings_wt), levels(rings_oe))))
  rings_wt <- factor(rings_wt, levels = all_levels)
  rings_oe <- factor(rings_oe, levels = all_levels)
  
  # Detection filtering (uses raw counts if available; otherwise data > 0)
  det_wt <- Matrix::rowMeans(X_wt > 0)
  det_oe <- Matrix::rowMeans(X_oe > 0)
  keep <- (det_wt >= min_detect_frac) & (det_oe >= min_detect_frac)
  
  genes <- rownames(X_wt)
  genes <- intersect(genes, rownames(X_oe))
  genes <- genes[keep[match(genes, rownames(X_wt))]]
  
  # Precompute ring index lists
  idx_wt <- split(seq_len(ncol(X_wt)), rings_wt)
  idx_oe <- split(seq_len(ncol(X_oe)), rings_oe)
  
  # Function to compute ring means for one gene across rings
  ring_means_gene <- function(X, idx_list, g) {
    v <- X[g, , drop = TRUE]
    sapply(idx_list, function(ii) if (length(ii) == 0) NA_real_ else mean(v[ii]))
  }
  
  out <- lapply(genes, function(g) {
    mw <- ring_means_gene(X_wt, idx_wt, g)
    mo <- ring_means_gene(X_oe, idx_oe, g)
    
    # Handle missing rings (no spots): set to 0 mass (after eps)
    mw[is.na(mw)] <- 0
    mo[is.na(mo)] <- 0
    
    pw <- profile_to_prob(mw, eps = eps)
    po <- profile_to_prob(mo, eps = eps)
    data.frame(
      gene = g,
      jsd = jsd(pw, po),
      mean_wt = mean(mw),
      mean_oe = mean(mo),
      det_wt = det_wt[match(g, rownames(X_wt))],
      det_oe = det_oe[match(g, rownames(X_oe))],
      stringsAsFactors = FALSE
    )
  })
  
  res <- do.call(rbind, out)
  res$log2fc_proxy <- res$mean_oe - res$mean_wt  # on log scale; proxy only
  res[order(res$jsd, decreasing = TRUE), ]
}

res_jsd <- ring_jsd_seurat(seu_wt, seu_oe, ring_col = "spa")
head(res_jsd, 20)


#######Plot Rings ##########
plot_ring_profiles <- function(seu_wt, seu_oe, gene,
                                  assay = DefaultAssay(seu_wt),
                                  layer = "data",
                                  ring_col = "ring") {
  
  # checks
  if (!ring_col %in% colnames(seu_wt@meta.data)) stop("WT: ring_col not in meta.data")
  if (!ring_col %in% colnames(seu_oe@meta.data)) stop("OE: ring_col not in meta.data")
  
  X_wt <- GetAssayData(seu_wt, assay = assay, layer = layer)
  X_oe <- GetAssayData(seu_oe, assay = assay, layer = layer)
  
  if (!gene %in% rownames(X_wt) || !gene %in% rownames(X_oe)) {
    warning("Gene not found in one of the objects: ", gene)
    return(invisible(NULL))
  }
  
  rings_wt <- seu_wt@meta.data[[ring_col]]
  rings_oe <- seu_oe@meta.data[[ring_col]]
  
  # drop NA rings (common source of split issues)
  keep_wt <- !is.na(rings_wt)
  keep_oe <- !is.na(rings_oe)
  
  if (sum(keep_wt) == 0) stop("WT: all ring labels are NA")
  if (sum(keep_oe) == 0) stop("OE: all ring labels are NA")
  
  rings_wt <- factor(rings_wt[keep_wt])
  rings_oe <- factor(rings_oe[keep_oe])
  
  # align ring levels
  all_levels <- sort(unique(c(levels(rings_wt), levels(rings_oe))))
  rings_wt <- factor(rings_wt, levels = all_levels)
  rings_oe <- factor(rings_oe, levels = all_levels)
  
  v_wt <- X_wt[gene, keep_wt, drop = TRUE]
  v_oe <- X_oe[gene, keep_oe, drop = TRUE]
  
  idx_wt <- split(seq_along(v_wt), rings_wt)
  idx_oe <- split(seq_along(v_oe), rings_oe)
  
  mw <- sapply(all_levels, function(lv) {
    ii <- idx_wt[[lv]]
    if (is.null(ii) || length(ii) == 0) NA_real_ else mean(v_wt[ii])
  })
  mo <- sapply(all_levels, function(lv) {
    ii <- idx_oe[[lv]]
    if (is.null(ii) || length(ii) == 0) NA_real_ else mean(v_oe[ii])
  })
  
  x <- seq_along(all_levels)
  y_range <- range(c(mw, mo), finite = TRUE)
  if (!is.finite(y_range[1])) y_range <- c(0, 0)
  
  plot(x, mw, type = "b", pch = 16, xaxt = "n",
       xlab = "Ring", ylab = "Mean log-normalized expr",
       main = gene, ylim = y_range)
  axis(1, at = x, labels = all_levels)
  lines(x, mo, type = "b", pch = 16)
  legend("topright", legend = c("WT", "OE"), lty = 1, pch = 16)
}


top_genes <- res_jsd$gene[1:12]
par(mfrow = c(3,4))
for (g in top_genes) plot_ring_profiles(seu_wt, seu_oe, g, ring_col = "spa")
par(mfrow = c(1,1))

plot_ring_profiles(seu_wt, seu_oe, top_genes[1], ring_col = "spa")


