getBaseATACPanel <- function(
  proj,
  fn,
  embedding = "UMAP"
) {
  baseGenes <- c(
    "CD3D",
    "CD4",
    "MPO",
    "CD8A",
    "MS4A1",
    "NCAM1",
    "TBX21",
    "EOMES",
    "SELL",
    "FAS",
    "FOXP3",
    "CD69",
    "CXCR5",
    "CCR5",
    "CD27",
    "CD28"
  )
  
  embed <- plotEmbedding(
    ArchRProj = proj,
    colorBy = "GeneScoreMatrix",
    name = baseGenes,
    embedding = embedding,
    baseSize = 5,
    imputeWeights = getImputeWeights(proj))
  
  embed2 <- lapply(seq_along(baseGenes), function(i) {
    g <- embed[[i]]
    g <- g + 
      labs(title = baseGenes[i]) + 
      theme(
        axis.title = element_blank(),
        plot.title = element_text(size = 8, margin = margin(7, 0, 0, 0)),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-10,5,0,-15),
        legend.box.spacing = unit(0, "pt"),
        legend.text = element_text(size = 5, margin = margin(t = -5)),
        legend.title = element_blank(),
        legend.key.size = unit(0.6, "lines"),
        legend.key.width = unit(1, "lines"),
        legend.position = "top",
        legend.justification = "right",
        plot.margin = unit(c(0,0,5,0), "pt"))
    
    return(g)
  })
  
  p <- patchwork::wrap_plots(embed2, ncol = 4)
  
  savePlot(plot = p, fn = fn, devices = c("png", "rds"), gheight = 10, gwidth = 8)
}


# initial cluster annotation
getBasePanelAndMarkers <- function(
  seu,
  tsa_catalog,
  plotFn,
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
  
  if (!("gene" %in% names(clusterMarkers))) {
    clusterMarkers$gene <- rownames(clusterMarkers)
  }
  
  clusterMarkersFilt <- clusterMarkers %>%
    left_join(tsa_catalog %>% dplyr::select(DNA_ID, cleanName), by = c("gene" = "DNA_ID")) %>% 
    filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    filter(avg_log2FC > 0) %>%
    mutate(piScore = avg_log2FC * -log10(p_val_adj)) %>%
    dplyr::arrange(desc(piScore), .by_group = TRUE)
  
  write.table(clusterMarkersFilt, file = tsvFn, row.names = FALSE, quote = FALSE, sep = "\t")
  
  # get base markers to help compare
  baseMarkers <- c(
    "A0034", #CD3
    "A0072", #CD4
    "A0046", #CD8
    "A0081", #CD14
    "A0083", #CD16
    "A0047", #CD56
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
    "A0088", #PD1
    "A0141", #CCR5
    "A0144", #CXCR5
    "A0149" #CD161
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
  
  
  savePlot(plot = baseRidgePlotGrid, fn = paste0(plotFn, "_separate"), devices = "rds")
  
  p <- do.call(cowplot::plot_grid, c(list(ncol = 5), baseRidgePlotGrid))
  savePlot(plot = p, fn = plotFn, devices = c("png", "rds"), gheight = 8, gwidth = 11)
}

# assign manual cluster annotation
assignManualAnnotation <- function(archrProj, annotFn, cluster, slotName = "manualClusterAnnot") {
  annotDf <- read.csv(annotFn, header = TRUE)
  
  manualAnnots <- annotDf$annotation
  names(manualAnnots) <- annotDf$cluster
  
  archrProj <- addCellColData(
    ArchRProj = archrProj,
    data = manualAnnots[getCellColData(archrProj, select = cluster)[, 1]],
    cells = archrProj$cellNames,
    name = "manualClusterAnnot",
    force = TRUE)
  
  return(archrProj)
}