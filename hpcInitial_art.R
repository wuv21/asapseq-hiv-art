suppressMessages(library(ArchR))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(rtracklayer))


addArchRThreads(threads = 4)
addArchRGenome("hg38")

# make vector of unwanted "chromosomes"/scaffolds
chrNames <- read.table("chrnames/all.chrnames.txt", header = FALSE)
chrNamesDiscard <- chrNames[grepl("(random|EBV|chrUn|chrM|chrHXB2)", chrNames[, 1]), 1]
chrNamesDiscard <- c(chrNamesDiscard, chrNames[grepl("^(KI|GL)", chrNames[, 1]), 1])
chrNamesDiscard <- c(chrNamesDiscard, chrNames[grepl("^(chrA08|chrA01)", chrNames[, 1]), 1])


# 10x file locations
samples_a08 <- c("A08_1", "A08_2", "A08_4")
fragmentsFN_a08 <- paste0("../20210620_asapseq_a08/count_out_a08Genome/", samples_a08, "/outs/fragments.tsv.gz")
# make arrow files
arrowfile_a08 <- createArrowFiles(inputFiles = fragmentsFN_a08,
  sampleNames = samples_a08,
  minTSS = 4,
  minFrags = 1000,
  addTileMat = TRUE,
  excludeChr = chrNamesDiscard,
  subThreading = TRUE, #some issue with hdf5
  force = FALSE, # if don't want to rewrite
  addGeneScoreMat = TRUE)

samples_a01 <- c("A01_1", "A01_2", "A01_3", "A01_4")
fragmentsFN_a01 <- paste0("../20210831_asapseq_a01/count_out_v2/", samples_a01, "/outs/fragments.tsv.gz")
arrowfile_a01 <- createArrowFiles(inputFiles = fragmentsFN_a01,
  sampleNames = samples_a01,
  minTSS = 4,
  minFrags = 1000,
  addTileMat = TRUE,
  excludeChr = chrNamesDiscard,
  subThreading = TRUE, #some issue with hdf5
  force = FALSE, # if don't want to rewrite
  addGeneScoreMat = TRUE)


combinedArrowFiles <- c(arrowfile_a08, arrowfile_a01)
initProjDir <- "A08A01_init"
qcFiltProjDir <- "A08A01_qcfilt"
if (!dir.exists(qcFiltProjDir)) {
  proj <- ArchRProject(
    ArrowFiles = combinedArrowFiles,
    outputDirectory = initProjDir,
    copyArrows = TRUE,
    showLogo = FALSE)

  idxPass <- BiocGenerics::which(proj$TSSEnrichment >= 9)
  cellsPass <- proj$cellNames[idxPass]
  projQcFilter <- proj[cellsPass, ]

  # filter using AMULET
  multCells <- lapply(c(samples_a08, samples_a01), function(x) {
    if (grepl("A08", x)) {
      amuletOutDir <- paste0("~/asapseq/20210620_asapseq_a08/amulet_out/", x)
    } else {
      amuletOutDir <- paste0("~/asapseq/20210831_asapseq_a01/amulet_out/", x)
    }

    amuletOut <- unlist(read.table(paste0(amuletOutDir, "/MultipletBarcodes_01.txt"), header = FALSE))
    amuletOut <- paste0(x, "#", amuletOut)

    return(amuletOut)
  })
  multCells <- unlist(multCells)

  singlets <- setdiff(projQcFilter$cellNames, multCells)
  projQcFilter <- projQcFilter[singlets, ]

  projQcFilter <- addIterativeLSI(
    ArchRProj = projQcFilter,
    useMatrix = "TileMatrix",
    name = "IterativeLSI",
    iterations = 2,
    clusterParams = list(resolution = c(0.5), n.start = 10),
    varFeatures = 25000,
    dimsToUse = 1:30
  )

  projQcFilter$individual <- ifelse(grepl("A01", projQcFilter$cellNames), "A01", "A08")
  projQcFilter <- addHarmony(
    ArchRProj = projQcFilter,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "individual"
  )

  projQcFilter <- saveArchRProject(
    ArchRProj = projQcFilter,
    outputDirectory = qcFiltProjDir,
    load = TRUE
  )
} else {
  projQcFilter <- loadArchRProject(
    path = qcFiltProjDir,
    showLogo = FALSE
  )
}

projQcFilter <- addClusters(
  input = projQcFilter,
  reducedDims = "Harmony",
  name = "Clusters",
  resolution = 0.3,
  force = TRUE
)

projQcFilter <- addUMAP(
  ArchRProj = projQcFilter,
  reducedDims = "Harmony",
  name = "UMAP",
  nNeighbors = 75,
  minDist = 0.1,
  metric = "cosine",
  force = TRUE
)

projQcFilter <- saveArchRProject(
  ArchRProj = projQcFilter,
  outputDirectory = qcFiltProjDir,
  load = TRUE
)

p1 <- plotEmbedding(ArchRProj = projQcFilter,
  colorBy = "cellColData",
  name = "Sample",
  embedding = "UMAP"
)

p2 <- plotEmbedding(ArchRProj = projQcFilter,
  colorBy = "cellColData",
  name = "Clusters",
  embedding = "UMAP"
)

plotPDF(p1, p2,
  name = "Plot-UMAP-Clusters.pdf",
  ArchRProj = projQcFilter,
  addDOC = FALSE,
  width = 7,
  height = 7
)

# markersGS <- getMarkerFeatures(
#   ArchRProj = projHighTSSLSI,
#   useMatrix = "GeneScoreMatrix",
#   groupBy = "Clusters",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   testMethod = "wilcoxon"
# )
# 
# saveRDS(markersGS, "rds/projTSSLSI_GS.rds")
# 
# 
