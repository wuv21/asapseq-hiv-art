#' Make Seurat antibody derived tags (ADT) object list
#'
#' @description
#' Saves and returns a RDS containing a list of Seurat objects by reading ADT files from kallisto bustools output.
#' 
#' @param samples A character vector of sample names that corresponds to the directory prefix
#' of the kallisto bustools output. For example, c("A01_1") will look for a directory that is
#' "A01_1_count_out", which uses the default "_count_out" suffix.
#' 
#' @param adtDir A character denoting where to ADT output directions are stored.
#' 
#' @param filename A character denoting where to save the RDS containing the RDS.
#' 
#' @param sampleSuffix A character value indicating the suffix of the directory containing the
#' ADT data. Default: _count_out.
#' 
#' @return A list of Seurat objects containing background filtered, normalized, and scaled ADT data.


createSeuList <- function(samples, adtDir, filename, sampleSuffix = "_count_out") {
  seuList <- lapply(samples, function(x) {
    bus_output <- read.table(file = glue("{adtDir}/{x}{sampleSuffix}/final.mtx"), header = FALSE, skip = 4)
    feats <- read.table(file = glue("{adtDir}/{x}{sampleSuffix}/final.genes.txt"), header = FALSE)
    cbcs <- read.table(file = glue("{adtDir}/{x}{sampleSuffix}/final.barcodes.revcomp.txt"), header = FALSE)
    
    adt_matrix <- Matrix::sparseMatrix(i = bus_output$V2,
      j = bus_output$V1,
      x = bus_output$V3)
    
    rownames(adt_matrix) <- feats$V1
    colnames(adt_matrix) <- cbcs$V1
    
    empty_calc <- emptyDrops(adt_matrix, lower = 500)
    adt_filt <- adt_matrix[, which(empty_calc$FDR < 0.01)]
    
    adt <- CreateSeuratObject(counts = adt_filt, assay = "tsa", project = x)
    
    adt <- NormalizeData(adt, normalization.method = "CLR", scale.factor = 10000, margin = 2)
    adt <- ScaleData(adt)
    
    DefaultAssay(adt) <- "tsa"
    
    return(adt)
  })
  
  saveRDS(seuList, file = filename)

  return(seuList)
}