suppressMessages(library(ArchR))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(rtracklayer))


addArchRThreads(threads = 4)
addArchRGenome("hg38")

# make vector of unwanted "chromosomes"/scaffolds
chrNames <- read.table("chrnames/all.chrnames.txt", header = FALSE)
chrNamesDiscard <- chrNames[grepl("(random|EBV|chrUn|chrM|chrHXB2)", chrNames[, 1]), 1]
chrNamesDiscard <- c(chrNamesDiscard, chrNames[grepl("^(KI|GL)", chrNames[, 1]), 1])

# 10x file locations
samples_c01 <- c("C01_1", "C01_2", "C01_3")
fragmentsFN_c01 <- paste0("../20211012_asapseq_c01/count_out_hxb2/", samples_c01, "/outs/fragments.tsv.gz")
arrowfile_c01 <- createArrowFiles(inputFiles = fragmentsFN_c01,
  sampleNames = samples_c01,
  minTSS = 4,
  minFrags = 1000,
  addTileMat = TRUE,
  excludeChr = chrNamesDiscard,
  subThreading = TRUE,
  force = FALSE,
  addGeneScoreMat = TRUE)

samples_c02 <- c("C02_1", "C02_2", "C02_3")
fragmentsFN_c02 <- paste0("../20211101_asapseq_c02/count_out_hxb2/", samples_c02, "/outs/fragments.tsv.gz")
arrowfile_c02 <- createArrowFiles(inputFiles = fragmentsFN_c02,
  sampleNames = samples_c02,
  minTSS = 4,
  minFrags = 1000,
  addTileMat = TRUE,
  excludeChr = chrNamesDiscard,
  subThreading = TRUE,
  force = FALSE,
  addGeneScoreMat = TRUE)


combinedArrowFiles <- c(arrowfile_c01, arrowfile_c02)
initProjDir <- "C01C02_init"
qcFiltProjDir <- "C01C02_qcfiltTSS6"
if (!dir.exists(qcFiltProjDir)) {
  proj <- ArchRProject(
    ArrowFiles = combinedArrowFiles,
    outputDirectory = initProjDir,
    copyArrows = TRUE,
    showLogo = FALSE)
  

  proj$individual <- ifelse(grepl("C01", proj$cellNames), "C01", "C02")
  
  #c01IdxPass <- BiocGenerics::which(proj$individual == "C01" & proj$TSSEnrichment >= 4)
  #c02IdxPass <- BiocGenerics::which(proj$individual == "C02" & proj$TSSEnrichment >= 8)

  idxPass <- BiocGenerics::which(proj$TSSEnrichment >= 6)
  cellsPass <- proj$cellNames[idxPass]
  projQcFilter <- proj[cellsPass, ]

  # filter using AMULET
  multCells <- lapply(c(samples_c01, samples_c02), function(x) {
    if (grepl("C01", x)) {
      amuletOutDir <- paste0("~/asapseq/20211012_asapseq_c01/amulet_out/", x)
    } else {
      amuletOutDir <- paste0("~/asapseq/20211101_asapseq_c02/amulet_out/", x)
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
