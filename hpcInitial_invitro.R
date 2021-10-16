suppressMessages(library(ArchR))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(rtracklayer))


addArchRThreads(threads = 4)
addArchRGenome("hg38")

# make vector of unwanted "chromosomes"/scaffolds
chrNames <- read.table("chrnames/all.chrnames.txt", header = FALSE)
chrNamesDiscard <- chrNames[grepl("(random|EBV|chrUn|chrM|chrHIV_suma)", chrNames[, 1]), 1]
chrNamesDiscard <- c(chrNamesDiscard, chrNames[grepl("^(KI|GL)", chrNames[, 1]), 1])
# chrNamesDiscard <- c(chrNamesDiscard, chrNames[grepl("^(chrA08|chrA01)", chrNames[, 1]), 1])


# 10x file locations
samples_nd497B <- c("ND497_B")
fragmentsFN_nd497B <- "../20210418_asapseq_pilot/cr2_ND497_B/outs/fragments.tsv.gz"
arrowfile_nd497B <- createArrowFiles(inputFiles = fragmentsFN_nd497B,
  sampleNames = samples_nd497B,
  minTSS = 4,
  minFrags = 1000,
  addTileMat = TRUE,
  excludeChr = chrNamesDiscard,
  subThreading = TRUE, #some issue with hdf5
  force = FALSE, # if don't want to rewrite
  addGeneScoreMat = TRUE)


initProjDir <- "ND497_B_init"
qcFiltProjDir <- "ND497_B_qcfilt"
if (!dir.exists(qcFiltProjDir)) {
  proj <- ArchRProject(
    ArrowFiles = arrowfile_nd497B,
    outputDirectory = initProjDir,
    copyArrows = TRUE,
    showLogo = FALSE)

  idxPass <- BiocGenerics::which(proj$TSSEnrichment >= 9)
  cellsPass <- proj$cellNames[idxPass]
  projQcFilter <- proj[cellsPass, ]

  # filter using AMULET
  multCells <- lapply(c(samples_nd497B), function(x) {
    amuletOutDir <- paste0("~/asapseq/20210418_asapseq_pilot/amulet_out/", x)

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

  projQcFilter$individual <- "ND497"

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
  reducedDims = "IterativeLSI",
  name = "Clusters",
  resolution = 0.8,
  force = TRUE
)

projQcFilter <- addUMAP(
  ArchRProj = projQcFilter,
  reducedDims = "IterativeLSI",
  name = "UMAP",
  nNeighbors = 30,
  minDist = 0.5,
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
