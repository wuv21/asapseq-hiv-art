# initial cluster annotation
findHIVDifferentialMarkers <- function(
  seu,
  tsa_catalog,
  identPos = NULL,
  identNeg = NULL,
  # cellsPos = NULL,
  # cellsNeg = NULL,
  featuresToUse = tsa_catalog[!tsa_catalog$isCtrl, "DNA_ID"],
  ident = "haystackOut",
  assay = "tsa",
  findMarkerMethod = "DESeq2",
  statusPos = "HIV+",
  statusNeg = "HIV-",
  ...
) {
  
  # Idents(seu) <- ident
  markers <- FindMarkers(
    object = seu,
    ident.1 = identPos,
    ident.2 = identNeg,
    assay = "tsa",
    test.use = findMarkerMethod,
    features = featuresToUse,
    ...)
  
  if (!"gene" %in% colnames(markers)) {
    markers$gene <- rownames(markers)
  }
  
  print(head(markers))
  
  filteredMarkers <- markers %>%
    left_join(tsa_catalog %>% dplyr::select(DNA_ID, cleanName), by = c("gene" = "DNA_ID")) %>% 
    mutate(p_val_adj = ifelse(p_val_adj == 0, .Machine$double.xmin, p_val_adj)) %>%
    filter(p_val_adj < 0.05) %>%
    # filter(avg_log2FC > 0) %>% 
    mutate(piScore = avg_log2FC * -log10(p_val_adj)) %>%
    arrange(desc(piScore), .by_group = TRUE) %>%
    mutate(Status = ifelse(avg_log2FC > 0, statusPos, statusNeg))
  
  # clean_VlnPlot(VlnPlot(seuMerge_matched, features = filteredHivPosMarkers$gene, combine = FALSE),
  #   newTitle = filteredHivPosMarkers$cleanName,
  #   finalplotCol = 4)
  
  return(filteredMarkers)
}






