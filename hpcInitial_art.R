suppressMessages(library(ArchR))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(rtracklayer))
suppressMessages(library(stringr))

set.seed(21)

addArchRThreads(threads = 2)
addArchRGenome("hg38")

# make vector of unwanted "chromosomes"/scaffolds
chrNames <- read.table("chrnames/all.chrnames.txt", header = FALSE)
chrNamesDiscard <- chrNames[grepl("(random|EBV|chrUn|chrM|chrHXB2)", chrNames[, 1]), 1]
chrNamesDiscard <- c(chrNamesDiscard, chrNames[grepl("^(KI|GL)", chrNames[, 1]), 1])
chrNamesDiscard <- c(chrNamesDiscard, chrNames[grepl("^(chrA08|chrA01|chrBEAT2045|chrA09)", chrNames[, 1]), 1])


# 10x file locations
samples_a08 <- c("A08_1", "A08_2", "A08_4")
fragmentsFN_a08 <- paste0("../asapseq-hiv-art-atacData/", samples_a08, "/fragments.tsv.gz")
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
fragmentsFN_a01 <- paste0("../asapseq-hiv-art-atacData/", samples_a01, "/fragments.tsv.gz")
arrowfile_a01 <- createArrowFiles(inputFiles = fragmentsFN_a01,
  sampleNames = samples_a01,
  minTSS = 4,
  minFrags = 1000,
  addTileMat = TRUE,
  excludeChr = chrNamesDiscard,
  subThreading = TRUE, #some issue with hdf5
  force = FALSE, # if don't want to rewrite
  addGeneScoreMat = TRUE)


samples_b45 <- c("BEAT045_A2", "BEAT045_A3", "BEAT045_B1", "BEAT045_B2")
sampleNames_b45 <- gsub("BEAT045", "B45", samples_b45)
fragmentsFN_b45 <- paste0("../asapseq-hiv-art-atacData/", samples_b45, "/fragments.tsv.gz")
arrowfile_b45 <- createArrowFiles(inputFiles = fragmentsFN_b45,
  sampleNames = sampleNames_b45,
  minTSS = 4,
  minFrags = 1000,
  addTileMat = TRUE,
  excludeChr = chrNamesDiscard,
  subThreading = TRUE, #some issue with hdf5
  force = FALSE, # if don't want to rewrite
  addGeneScoreMat = TRUE)

samples_a09 <- c("A09_pre_1", "A09_pre_2", "A09_pre_3", "A09_pre_4", "A09_pre_5", "A09_pre_6",
                "A09_post_2", "A09_post_3", "A09_post_4", "A09_post_5", "A09_post_6")
fragmentsFN_a09 <- paste0("../asapseq-hiv-art-atacData/", samples_a09, "/fragments.tsv.gz")
arrowfile_a09 <- createArrowFiles(inputFiles = fragmentsFN_a09,
                                  sampleNames = samples_a09,
                                  minTSS = 4,
                                  minFrags = 1000,
                                  addTileMat = TRUE,
                                  excludeChr = chrNamesDiscard,
                                  subThreading = TRUE, #some issue with hdf5
                                  force = FALSE, # if don't want to rewrite
                                  addGeneScoreMat = TRUE)



combinedArrowFiles <- c(arrowfile_a08, arrowfile_a01, arrowfile_b45, arrowfile_a09)
initProjDir <- "A08A01B45A09_init"
qcFiltProjDir <- "A08A01B45A09_qcfiltTSS8"
if (!dir.exists(qcFiltProjDir)) {
  proj <- ArchRProject(
    ArrowFiles = combinedArrowFiles,
    outputDirectory = initProjDir,
    copyArrows = TRUE,
    showLogo = FALSE)

  idxPass <- BiocGenerics::which(proj$TSSEnrichment >= 8)
  cellsPass <- proj$cellNames[idxPass]
  projQcFilter <- proj[cellsPass, ]

  # filter using AMULET
  multCells <- lapply(c(samples_a08, samples_a01, samples_b45, samples_a09), function(x) {

    amuletOutDir <- paste0("../asapseq-hiv-art-atacData/", x)

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

  projQcFilter$individual <- stringr::str_match(projQcFilter$cellNames, "\\w\\d{2}(_pre|_post)*")[, 1]

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
  dropCells = TRUE,
  outputDirectory = qcFiltProjDir,
  load = TRUE
)
