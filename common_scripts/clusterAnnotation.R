# initial cluster annotation
getBasePanelAndMarkers <- function(seu,
  tsa_catalog,
  pngFn,
  tsvFn,
  featuresToUse = rownames(tsa_catalog[!tsa_catalog$isCtrl, "DNA_ID"]),
  findMarkerMethod = "wilcox",
  clusterName = "atacClusters") {
  
  Idents(seu) <- clusterName
  clusterMarkers <- FindAllMarkers(seu,
    assay = "tsa",
    slot = "data",
    test.use = findMarkerMethod,
    features = featuresToUse)
  
  clusterMarkersFilt <- clusterMarkers %>%
    left_join(tsa_catalog %>% dplyr::select(DNA_ID, cleanName), by = c("gene" = "DNA_ID")) %>% 
    filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    filter(avg_log2FC > 0) %>%
    mutate(piScore = avg_log2FC * -log10(p_val_adj)) %>%
    arrange(desc(piScore), .by_group = TRUE)
  
  write.table(clusterMarkersFilt, file = tsvFn, row.names = FALSE, quote = FALSE, sep = "\t")
  
  # get base markers to help compare
  baseMarkers <- c(
    "A0034", #CD3
    "A0072", #CD4
    "A0046", #CD8
    "A0081", #CD14
    "A0083", #CD16
    "A0050", #CD19
    "A0087", #CD45RO
    "A0063", #CD45RA
    "A0154", #CD27
    "A0386", #CD28
    "A0156", #CD95
    "A0390", #CD127
    "A0146", #CD69
    "A0085", #CD25
    "A0089", #TIGIT
    "A0088" #PD1
  )
  
  baseRidgePlot <- RidgePlot(seu,
    features = baseMarkers,
    group.by = clusterName,
    combine = FALSE)
  
  baseRidgePlotGrid <- lapply(baseRidgePlot, function(x) {
    tsaID <- x$labels$title
    tsaCDName <- tsa_catalog[tsa_catalog$DNA_ID == tsaID, "cleanName"]
    
    x +
      guides(color = "none", fill = "none") +
      labs(title = glue("{tsaCDName} ({tsaID})")) +
      theme(axis.title = element_blank(),
        axis.text = element_text(size = 6),
        panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(size = 8))
  })
  
  cowplot::save_plot(filename = pngFn,
    plot = do.call(cowplot::plot_grid, c(list(ncol = 5), baseRidgePlotGrid)),
    base_height = 7,
    base_width = 10,
    bg = "#ffffff",
    dpi = "retina")
}

# assign manual cluster annotation
assignManualAnnotation <- function(archrProj, annotFn, slotName = "manualClusterAnnot") {
  annotDf <- read.csv(annotFn, header = TRUE)
  
  manualAnnots <- annotDf$annotation
  names(manualAnnots) <- annotDf$cluster
  
  archrProj <- addCellColData(
    ArchRProj = archrProj,
    data = manualAnnots[archrProj$Clusters],
    cells = archrProj$cellNames,
    name = "manualClusterAnnot",
    force = TRUE)
  
  return(archrProj)
}