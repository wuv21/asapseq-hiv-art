.testMarkerSCCustom <- function(
  ArrowFiles = NULL,
  matchObj = NULL,
  group = NULL,
  testMethod = "ttest",
  useMatrix = NULL,
  threads = 1,
  featureDF,
  binarize = FALSE,
  normFactors = NULL,
  logFile = NULL
){
  
  matchx <- matchObj[[1]][[group]]
  cellsx <- matchObj[[2]]$cells[matchx$cells]
  bgdx <- matchObj[[2]]$cells[matchx$bgd]
  
  if(!is.null(normFactors)){
    cellNF <- as.numeric(normFactors[cellsx,1])
    bgdNF <- as.numeric(normFactors[bgdx,1])
  }
  
  #Add RowNames for Check at the end
  rownames(featureDF) <- paste0("f", seq_len(nrow(featureDF)))
  seqnames <- unique(featureDF$seqnames)
  
  .logThis(cellsx, paste0(group, "_cellsx"), logFile = logFile)
  .logThis(bgdx, paste0(group, "_bgdx"), logFile = logFile)
  
  pairwiseDF <- lapply(seq_along(seqnames), function(y){
    
    .logMessage(sprintf("Pairwise Test %s : Seqnames %s", group, seqnames[y]), logFile = logFile)
    featureDFy <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% seqnames[y]), ]
    
    if(length(c(cellsx, bgdx)) == 0){
      stop(paste0("Cells in foreground and background are 0 for group = ", group))
    }
    
    scMaty <- suppressMessages(.getPartialMatrix(
      ArrowFiles, 
      featureDF = featureDFy, 
      threads = threads, 
      useMatrix = useMatrix,
      cellNames = c(cellsx, bgdx),
      progress = FALSE
    ))
    scMaty <- .checkSparseMatrix(scMaty, length(c(cellsx, bgdx)))
    
    .logThis(scMaty, paste0(group, "_", seqnames[y], "_scMaty"), logFile = logFile)
    rownames(scMaty) <- rownames(featureDFy)
    
    if(binarize){
      scMaty@x[scMaty@x > 0] <- 1
    }
    
    args <- list()
    
    if(!is.null(normFactors)){
      args$mat1 <- Matrix::t(Matrix::t(scMaty[, cellsx, drop = FALSE]) * cellNF)
      args$mat2 <- Matrix::t(Matrix::t(scMaty[, bgdx, drop = FALSE]) * bgdNF)
    }else{
      args$mat1 <- scMaty[, cellsx, drop = FALSE]
      args$mat2 <- scMaty[, bgdx, drop = FALSE]
    }
    
    if(tolower(testMethod) == "wilcoxon"){
      
      o <- tryCatch({
        .suppressAll(do.call(.sparseMatWilcoxon, args))
      }, error = function(e){
        errorList <- args
        .logError(e, fn = ".sparseMatWilcoxon", info = seqnames[y], errorList = errorList, logFile = logFile)
      })
      
    }else if(tolower(testMethod) == "ttest"){
      
      o <- tryCatch({
        .suppressAll(do.call(.sparseMatTTest, args))
      }, error = function(e){
        errorList <- args
        .logError(e, fn = ".sparseMatTTest", info = seqnames[y], errorList = errorList, logFile = logFile)
      })
      
    }else if(tolower(testMethod) == "binomial"){
      
      if(!is.null(normFactors)){
        .logStop("Normfactors cannot be used with a binomial test!", logFile = logFile)
      }
      
      if(!binarize){
        .logStop("Binomial test requires binarization!", logFile = logFile)
      }
      
      o <- tryCatch({
        .suppressAll(do.call(.sparseMatBinomTest, args))
      }, error = function(e){
        errorList <- args
        .logError(e, fn = ".sparseMatBinomTest", info = seqnames[y], errorList = errorList, logFile = logFile)
      })
      
    }else{
      
      .logStop("Error Unrecognized Method!", logFile = logFile)
      
    }
    
    .logThis(o, paste0(group, "_", seqnames[y], "_diffResult"), logFile = logFile)
    
    o
    
  }) %>% Reduce("rbind", .)
  
  idxFilter <- rowSums(pairwiseDF[, c("mean1", "mean2")]) != 0 & rowSums(is.na(pairwiseDF[, c("mean1", "mean2")])) == 0 # updated this line.
  pairwiseDF$fdr <- NA
  pairwiseDF$fdr[idxFilter] <- p.adjust(pairwiseDF$pval[idxFilter], method = "fdr")
  pairwiseDF <- pairwiseDF[rownames(featureDF), , drop = FALSE]
  pairwiseDF
  
}