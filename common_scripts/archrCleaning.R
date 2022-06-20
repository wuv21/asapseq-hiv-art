######
# Filter small clusters based on proportion
######
filterClusterByProp <- function(proj, proportion, cluster) {
  # filter cells that are in small clusters
  clusters <- getCellColData(proj, select = cluster)[, 1]
  
  clusterSize <- table(clusters) / length(proj$cellNames) * 100
  clusterSizeFilt <- clusterSize >= proportion
  
  print(paste0("filtering out outliers: ", length(proj$cellNames) - sum(clusterSizeFilt)))
  
  idxClusterPass <- BiocGenerics::which(clusterSizeFilt[clusters])
  cellsClusterPass <- proj$cellNames[idxClusterPass]
  proj_clustFilt <- proj[cellsClusterPass, ]
  
  return(proj_clustFilt)
}

addHaystackData <- function(proj, haystackParentDir, haystackSamples) {
  haystackList <- lapply(haystackSamples, function(i) {
    df <- read.csv(glue("{haystackParentDir}/{i}/viralFrags.tsv"), header = TRUE, sep = "\t")
    df <- df %>%
      filter(alreadyRecordedInIntegration != "True")
    
    if (nrow(df) == 0) {
      return(NULL)
    }
    
    intSiteFragsFn <- glue("{haystackParentDir}/{i}/integrationSites_viralFrags.tsv")
    if (file.exists(intSiteFragsFn)) {
      df2 <- read.csv(intSiteFragsFn, header = TRUE, sep = "\t")
      
      if (nrow(df2) > 0) {
        df <- bind_rows(df, df2)
      }

    }
    return(df)
  })
  names(haystackList) <- names(haystackSamples)
  
  haystackDf <- bind_rows(haystackList, .id = "sample")
  haystackDf <- haystackDf %>%
    mutate(newCbc = paste0(sample, "#", cbc))
  
  haystackPosCells <- haystackDf %>%
    group_by(sample, newCbc) %>%
    dplyr::summarize(totalFrags = n())
  
  matchedHaystackArchRCells <- intersect(haystackPosCells$newCbc, proj$cellNames)
  
  haystackPosCellsViralFrags <- haystackDf %>%
    filter(newCbc %in% matchedHaystackArchRCells)
  
  haystackOut <- rep(FALSE, length(proj$cellNames))
  names(haystackOut) <- proj$cellNames
  haystackOut[matchedHaystackArchRCells] <- TRUE
  
  proj$haystackOut <- haystackOut
  
  return(list(newProj = proj, allViralFrags = haystackDf, filteredviralFrags = haystackPosCellsViralFrags))
}

######
# Find ADT markers similar to background
######
returnSimilarToIsotype <- function(tsa_catalog, isoComparison, pThresh = 0.05, nSimilar = 4) {
  isoComparisonSum <- rowSums(isoComparison > pThresh)
  adtSimilar <- tsa_catalog[!tsa_catalog$isCtrl, ] %>%
    mutate(similarToControl = isoComparisonSum[DNA_ID]) %>%
    dplyr::select(DNA_ID, cleanDescriptionFull, similarToControl) %>%
    dplyr::filter(similarToControl > nSimilar)
  
  return(adtSimilar$DNA_ID)
}

######
# get ADT markers not similar to background
######
returnADTNonSimilar <- function(seu, tsa_catalog, isoSimilar) {
  return(rownames(seu)[!tsa_catalog$isCtrl & !(rownames(seu) %in% isoSimilar)])
}


######
# get intersecting cell barcodes
######
getIntersectBarcodes <- function(seu, archrproj) {
  intersectCbc <- intersect(Cells(seu), archrproj$cellNames)
  
  return(intersectCbc)
}
