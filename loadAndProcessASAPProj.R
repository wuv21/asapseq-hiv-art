.getMatchedProjName <- function(proj, suffix = "_matched") {
  return(paste0(getOutputDirectory(proj), suffix))
}

loadProcessedProj <- function(
  proj, # archr project
  adt, # seurat object containing adt data
  tsaCatalog, # data frame of tsa catalog
  generateNewUmapSettings = NULL, # list of parameters for addUMAP to run on archr project
  generateNewClusterSettings = NULL, # list of parameters for addClusters to run on archr proj
  umapName = ifelse(!is.null(generateNewUmapSettings), generateNewUmapSettings[["name"]], "UMAP"),
  clusterName = ifelse(!is.null(generateNewClusterSettings), generateNewClusterSettings[["name"]], "Clusters"),
  runPreAnnot = TRUE,
  runPreAnnotIsoComps = NULL,
  plotGraphs = TRUE, # plot graphs,
  prefixForGraphs = deparse(substitute(proj)), # name of files to save
  force = FALSE # remake processed folder
) {
  
  # check if processed folder exists
  projMatchedDir <- .getMatchedProjName(proj)
  if (force | !dir.exists(projMatchedDir)) {
    print("Did not find processed dir or force = TRUE. Beginning processing")
    
    # get intersecting cbc
    intersectCbc <- getIntersectBarcodes(adt, proj)
    proj_matched <- proj[intersectCbc, ]
    
    print("UMAP in process")
    # add umap if requested
    if (!is.null(generateNewUmapSettings)) {
      proj_matched <- do.call(addUMAP, c(list("ArchRProj" = proj_matched), generateNewUmapSettings))
    }
    
    print("Clustering in process")
    # add cluster if requested
    if (!is.null(generateNewClusterSettings)) {
      proj_matched <- do.call(addClusters, c(list("input" = proj_matched), generateNewClusterSettings))
    }
    
    # adjust cluster
    proj_matched <- filterClusterByProp(proj_matched, proportion = 0.5, cluster = clusterName)

  } else {
    proj_matched <- loadArchRProject(projMatchedDir, force = TRUE, showLogo = FALSE)
  }
  
  adt_matched <- subset(adt, cells = proj_matched$cellNames)
  
  if (plotGraphs) {
    print("Plotting in process")
    
    plotUmap(
      proj_matched,
      fn = paste0(prefixForGraphs, "_umap_unannotCluster"),
      colorBy = clusterName,
      embedding = umapName,
      colorScheme = scale_color_manual(values = multiHueColorPalette))
    
    plotUmap(
      proj_matched,
      fn = paste0(prefixForGraphs, "_umap_haystackout"),
      colorBy = "haystackOut",
      colorLabelBy = NULL,
      colorScheme = NULL,
      bringToTop = TRUE,
      embedding = umapName,
      propInLegend = TRUE)
    
    if (length(unique(proj_matched$individual)) > 1) {
      plotUmap(
        proj_matched,
        colorBy = "individual",
        embedding = umapName,
        devices = "png",
        fn = paste0(prefixForGraphs, "_umap_indivdiual"))
    }
  }
  
  isoSimilar <- returnSimilarToIsotype(tsaCatalog, runPreAnnotIsoComps, 0.05, 4)
  adtToUse <- returnADTNonSimilar(adt, tsaCatalog, isoSimilar)
  adt_matched$atacClusters <- getCellColData(proj_matched, clusterName)[, 1]
  
  if (runPreAnnot) {
    print("Pre annotation in process")
    
    proj_matched <- addImputeWeights(proj_matched)
    
    getBaseATACPanel(proj_matched,
      paste0(prefixForGraphs, "_imputePanel"),
      embedding = umapName)
    
    getBasePanelAndMarkers(
      seu = adt_matched,
      tsa_catalog = tsa_catalog,
      featuresToUse = adtToUse,
      findMarkerMethod = "wilcox",
      plotFn = paste0(prefixForGraphs, "_basepanel"), 
      tsvFn = paste0("outs/tsv/", prefixForGraphs, "_initClusterMarkers.tsv"))
  }
  
  proj_matched <- saveArchRProject(proj_matched, outputDirectory = projMatchedDir)
  
  return(list(
    "adtToUse" = adtToUse,
    "archrProj" = proj_matched,
    "adtSeu" = adt_matched
  ))
}