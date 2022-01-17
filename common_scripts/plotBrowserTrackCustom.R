plotBrowserTrack2 <- function(
  ArchRProj = NULL, 
  region = NULL, 
  groupBy = "Clusters",
  useGroups = NULL, 
  plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
  sizes = c(10, 1.5, 3, 4),
  features = getPeakSet(ArchRProj),
  loops = getCoAccessibility(ArchRProj),
  geneSymbol = NULL,
  useMatrix = NULL,
  log2Norm = TRUE,
  upstream = 50000,
  downstream = 50000,
  tileSize = 250, 
  minCells = 25,
  normMethod = "ReadsInTSS",
  threads = getArchRThreads(), 
  ylim = NULL,
  pal = NULL,
  baseSize = 7,
  scTileSize = 0.5,
  scCellsMax = 100,
  borderWidth = 0.4,
  tickWidth = 0.4,
  facetbaseSize = 7,
  geneAnnotation = getGeneAnnotation(ArchRProj),
  title = "",
  verbose = TRUE,
  logFile = createLogFile("plotBrowserTrack")
){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProj")
  .validInput(input = region, name = "region", valid = c("granges","null"))
  .validInput(input = groupBy, name = "groupBy", valid = "character")
  .validInput(input = useGroups, name = "useGroups", valid = c("character", "null"))
  .validInput(input = plotSummary, name = "plotSummary", valid = "character")
  .validInput(input = sizes, name = "sizes", valid = "numeric")
  .validInput(input = features, name = "features", valid = c("granges", "grangeslist", "null"))
  .validInput(input = loops, name = "loops", valid = c("granges", "grangeslist", "null"))
  .validInput(input = geneSymbol, name = "geneSymbol", valid = c("character", "null"))
  .validInput(input = useMatrix, name = "useMatrix", valid = c("character", "null"))
  .validInput(input = log2Norm, name = "log2Norm", valid = c("boolean"))
  .validInput(input = upstream, name = "upstream", valid = c("integer"))
  .validInput(input = downstream, name = "downstream", valid = c("integer"))
  .validInput(input = tileSize, name = "tileSize", valid = c("integer"))
  .validInput(input = minCells, name = "minCells", valid = c("integer"))
  .validInput(input = normMethod, name = "normMethod", valid = c("character"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = ylim, name = "ylim", valid = c("numeric", "null"))
  .validInput(input = pal, name = "pal", valid = c("palette", "null"))
  .validInput(input = baseSize, name = "baseSize", valid = "numeric")
  .validInput(input = scTileSize, name = "scTileSize", valid = "numeric")
  .validInput(input = scCellsMax, name = "scCellsMax", valid = "integer")
  .validInput(input = borderWidth, name = "borderWidth", valid = "numeric")
  .validInput(input = tickWidth, name = "tickWidth", valid = "numeric")
  .validInput(input = facetbaseSize, name = "facetbaseSize", valid = "numeric")
  geneAnnotation <- .validGeneAnnotation(geneAnnotation)
  .validInput(input = title, name = "title", valid = "character")
  
  tstart <- Sys.time()
  .startLogging(logFile=logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "plotBrowserTrack Input-Parameters", logFile = logFile)
  
  ##########################################################
  # Get Region Where Plot Will Occur (GenomicRanges)
  ##########################################################
  .logDiffTime("Validating Region", t1=tstart, verbose=verbose, logFile=logFile)
  if(is.null(region)){
    if(!is.null(geneSymbol)){
      region <- geneAnnotation$genes
      region <- region[which(tolower(mcols(region)$symbol) %in% tolower(geneSymbol))]
      region <- region[order(match(tolower(mcols(region)$symbol), tolower(geneSymbol)))]
      print(region)
      region <- resize(region, 1, "start")
      strand(region) <- "*"
      region <- extendGR(region, upstream = upstream, downstream = downstream)
    }
  }
  region <- .validGRanges(region)
  .logThis(region, "region", logFile = logFile)
  
  if(is.null(geneSymbol)){
    useMatrix <- NULL
  }
  
  if(!is.null(useMatrix)){
    featureMat <- .getMatrixValues(
      ArchRProj = ArchRProj,
      matrixName = useMatrix,
      name = mcols(region)$symbol
    )
    if(log2Norm){
      featureMat <- log2(featureMat + 1) 
    }
    featureMat <- data.frame(t(featureMat))
    featureMat$Group <- getCellColData(ArchRProj, groupBy, drop = FALSE)[rownames(featureMat), 1]
  }
  
  ggList <- lapply(seq_along(region), function(x){
    
    plotList <- list()
    
    ##########################################################
    # Bulk Tracks
    ##########################################################
    if("bulktrack" %in% tolower(plotSummary)){
      .logDiffTime(sprintf("Adding Bulk Tracks (%s of %s)",x,length(region)), t1=tstart, verbose=verbose, logFile=logFile)
      plotList$bulktrack <- .bulkTracks(
        ArchRProj = ArchRProj, 
        region = region[x], 
        tileSize = tileSize, 
        groupBy = groupBy,
        threads = threads, 
        minCells = minCells,
        pal = pal,
        ylim = ylim,
        baseSize = baseSize,
        borderWidth = borderWidth,
        tickWidth = tickWidth,
        facetbaseSize = facetbaseSize,
        normMethod = normMethod,
        geneAnnotation = geneAnnotation,
        title = title,
        useGroups = useGroups,
        tstart = tstart,
        logFile = logFile) + theme(plot.margin = unit(c(0.35, 0.75, 0.35, 0.75), "cm"))
    }
    
    ##########################################################
    # Bulk Tracks
    ##########################################################
    if("sctrack" %in% tolower(plotSummary)){
      .logDiffTime(sprintf("Adding SC Tracks (%s of %s)",x,length(region)), t1=tstart, verbose=verbose, logFile=logFile)
      plotList$sctrack <- .scTracks(
        ArchRProj = ArchRProj, 
        region = region[x], 
        tileSize = tileSize, 
        groupBy = groupBy,
        threads = threads, 
        minCells = 5,
        maxCells = scCellsMax,
        pal = pal,
        baseSize = baseSize,
        borderWidth = borderWidth,
        tickWidth = tickWidth,
        scTileSize = scTileSize,
        facetbaseSize = facetbaseSize,
        geneAnnotation = geneAnnotation,
        title = title,
        useGroups = useGroups,
        tstart = tstart,
        logFile = logFile) + theme(plot.margin = unit(c(0.35, 0.75, 0.35, 0.75), "cm"))
    }
    
    ##########################################################
    # Feature Tracks
    ##########################################################
    if("featuretrack" %in% tolower(plotSummary)){
      if(!is.null(features)){
        .logDiffTime(sprintf("Adding Feature Tracks (%s of %s)",x,length(region)), t1=tstart, verbose=verbose, logFile=logFile)
        plotList$featuretrack <- .featureTracks(
          features = features, 
          region = region[x], 
          facetbaseSize = facetbaseSize,
          hideX = TRUE, 
          title = "Peaks",
          logFile = logFile) + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
      }
    }
    
    ##########################################################
    # Feature Tracks
    ##########################################################
    if("looptrack" %in% tolower(plotSummary)){
      if(!is.null(loops)){
        .logDiffTime(sprintf("Adding Loop Tracks (%s of %s)",x,length(region)), t1=tstart, verbose=verbose, logFile=logFile)
        plotList$looptrack <- .loopTracks(
          loops = loops, 
          region = region[x], 
          facetbaseSize = facetbaseSize,
          hideX = TRUE, 
          hideY = TRUE,
          title = "Loops",
          logFile = logFile) + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
      }
    }
    
    ##########################################################
    # Gene Tracks
    ##########################################################
    if("genetrack" %in% tolower(plotSummary)){
      .logDiffTime(sprintf("Adding Gene Tracks (%s of %s)",x,length(region)), t1=tstart, verbose=verbose, logFile=logFile)
      plotList$genetrack <- .geneTracks(
        geneAnnotation = geneAnnotation, 
        region = region[x], 
        facetbaseSize = facetbaseSize,
        title = "Genes",
        logFile = logFile) + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
    }
    
    ##########################################################
    # Time to plot
    ##########################################################
    plotSummary <- tolower(plotSummary)
    names(sizes) <- plotSummary
    sizes <- sizes[order(plotSummary)]
    plotSummary <- plotSummary[order(plotSummary)]
    
    # nullSummary <- unlist(lapply(seq_along(plotSummary), function(x) is.null(eval(parse(text=paste0("plotList$", plotSummary[x]))))))
    # if(any(nullSummary)){
    #   sizes <- sizes[-which(nullSummary)]
    # }
    sizes <- sizes[tolower(names(plotList))]
    
    if(!is.null(useMatrix)){
      
      suppressWarnings(.combinedFeaturePlot(
        plotList = plotList,
        log2Norm = log2Norm,
        featureMat = featureMat,
        feature = region[x]$symbol[[1]],
        useMatrix = useMatrix,
        pal = pal,
        sizes = sizes,
        baseSize = baseSize,
        facetbaseSize = facetbaseSize,
        borderWidth = borderWidth,
        tickWidth = tickWidth
      ))
      
    }else{
      
      .logThis(names(plotList), sprintf("(%s of %s) names(plotList)",x,length(region)), logFile=logFile)
      .logThis(sizes, sprintf("(%s of %s) sizes",x,length(region)), logFile=logFile)
      #.logThis(nullSummary, sprintf("(%s of %s) nullSummary",x,length(region)), logFile=logFile)
      .logDiffTime("Plotting", t1=tstart, verbose=verbose, logFile=logFile)
      
      tryCatch({
        # suppressWarnings(ggAlignPlots(plotList = plotList, sizes=sizes, draw = FALSE))
        return(plotList)
      }, error = function(e){
        .logMessage("Error with plotting, diagnosing each element", verbose = TRUE, logFile = logFile)
        for(i in seq_along(plotList)){
          tryCatch({
            print(plotList[[i]])
          }, error = function(f){
            .logError(f, fn = names(plotList)[i], info = "", errorList = NULL, logFile = logFile)
          })
        }
        .logError(e, fn = "ggAlignPlots", info = "", errorList = NULL, logFile = logFile)
      })
      
    }
    
  })
  
  if(!is.null(mcols(region)$symbol)){
    names(ggList) <- mcols(region)$symbol
  }else{
    if(length(ggList) == 1){
      ggList <- ggList[[1]]
    }
  }
  
  .endLogging(logFile=logFile)
  
  ggList
  
}