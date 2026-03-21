SCT_splitwise <- function(data, split.by = "group", assay ="Spatial",
                          slot = "counts"){
  data.list <- SplitObject(data, split.by = split.by)
  data.norm.list <- list()
  dds.list <- list()
  var_feat.list <- list()

  for (sample in data.list){
    sample <- SCTransform(sample, assay = assay, new.assay.name = "SCT",
                          return.only.var.genes = FALSE, variable.features.n = 3000)
    #vars.to.regress = "Wnt1") # If Wnt1 was orders of magnitudes higher
    data.norm.list <- append(data.norm.list, sample)
    var_feat.sample <- list(VariableFeatures(sample, assay = "SCT"))
    var_feat.list <- c(var_feat.list, var_feat.sample)
  }

  # Proper intersection across all groups
  var_feat <- Reduce(intersect, var_feat.list)

  # Merge ALL groups (not only the first two)
  merged <- Reduce(function(x, y) merge(x, y), data.norm.list)

  DefaultAssay(merged) <- "SCT"
  VariableFeatures(merged) <- var_feat
  merged <- PrepSCTFindMarkers(merged)

  # Keep original image (minimal change from your approach)
  merged@images$slice1 <- data@images$slice1
  merged@images$slice1.2 <- NULL

  return(merged)
}
