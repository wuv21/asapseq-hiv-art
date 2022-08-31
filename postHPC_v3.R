###############################################################################
# pkg and custom scripts
###############################################################################
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ArchR))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(rtracklayer))
library(patchwork)
library(Seurat)
library(glue)
library(RColorBrewer)
library(ggrastr)
library(ggpubr)
library(biomaRt)
library(ComplexHeatmap)
library(BSgenome.Hsapiens.UCSC.hg38)

source("common_scripts/graphs.R")
source("common_scripts/archrCleaning.R")
source("common_scripts/clusterAnnotation.R")
source("common_scripts/differentialExpression.R")
source("common_scripts/exportTsv.R")
source("common_scripts/colorPalettes.R")

source("loadAndProcessASAPProj.R")

#https://stackoverflow.com/questions/23279904/modifying-an-r-package-function-for-current-r-session-assigninnamespace-not-beh
source("common_scripts/plotBrowserTrackCustom.R")
source("common_scripts/testSCMarkerCustom.R")

select <- dplyr::select
slice <- dplyr::slice

tmpFun <- get("plotBrowserTrack", envir = asNamespace("ArchR"))
environment(plotBrowserTrack2) <- environment(tmpFun)
attributes(plotBrowserTrack2) <- attributes(tmpFun)

environment(.testMarkerSCCustom) <- asNamespace('ArchR')
assignInNamespace(".testMarkerSC", .testMarkerSCCustom, ns = "ArchR")

###############################################################################
# analysis settings
###############################################################################
addArchRThreads(threads = 3) # can adjust if more RAM is available (higher threads require more)
set.seed(21) # for reproducibility

###############################################################################
# load ADT files
###############################################################################
tsa_catalog <- readRDS("rds/tsa_catalog.rds")
load("rds/invitro_seuratMerged.RData")
load("rds/ART_seuratMerged.RData")

###############################################################################
# load in iv project from preHPC 
###############################################################################
projInVitro <- loadArchRProject(path = "ND497_B_qcfiltTSS8/", showLogo = FALSE)

haystackInVitro <- addHaystackData(proj = projInVitro,
  haystackParentDir = "data/hiv-haystack",
  haystackSamples = c("ND497_B" = "ND497_B_v3"))

projInVitro <- haystackInVitro$newProj


###############################################################################
# load in art project from preHPC 
###############################################################################
projART <- loadArchRProject(path = "A08A01B45A09_qcfiltTSS8/", showLogo = FALSE)

# add haystack data
# need to edit because samples is edited
uniqSamples <- unique(projART$Sample)
haystackSamplePairingART <- sapply(uniqSamples, function(i) {
  if (grepl("A01_post", i)) {
    newSampleName <- gsub("_post", "", i)
    newSampleName <- paste0(newSampleName, "_v3WithEnv_filtered")
  } else if (grepl("A08_post", i)) {
    newSampleName <- gsub("_post", "", i)
    newSampleName <- paste0(newSampleName, "_v3")
  } else {
    newSampleName <- gsub("B45", "BEAT045", i)
    newSampleName <- paste0(newSampleName, "_consensus")
  }
  
  return(newSampleName)
})

haystackART <- addHaystackData(proj = projART,
  haystackParentDir = "data/hiv-haystack",
  haystackSamples = haystackSamplePairingART)

projART <- haystackART$newProj

###############################################################################
# iv processing
###############################################################################

ivPreProcessed <- loadProcessedProj(
  proj = projInVitro,
  adt = adtInVitro,
  tsaCatalog = tsa_catalog,
  generateNewUmapSettings = list(
    "reducedDims" = "IterativeLSI",
    "name" = "UMAP2",
    "nNeighbors" = 300,
    "force" = TRUE,
    "minDist" = 0.4),
  generateNewClusterSettings = list(
    "method" = "Seurat",
    "name" = "Clusters2",
    "resolution" = 0.8
  ),
  runPreAnnot = TRUE,
  runPreAnnotIsoComps = isoComparisonsInvitro,
  prefixForGraphs = "invitro"
)

projInvitro_matched <- ivPreProcessed$archrProj
adtInvitro_matched <- ivPreProcessed$adtSeu
invitroAdtToUse <- ivPreProcessed$adtToUse

rm(ivPreProcessed)
gc()

###############################################################################
# art processing
###############################################################################

artPreProcessed <- loadProcessedProj(
  proj = projART,
  adt = adtART,
  tsaCatalog = tsa_catalog,
  generateNewUmapSettings = list(
    "reducedDims" = "Harmony",
    "name" = "UMAP2",
    "nNeighbors" = 100,
    "force" = TRUE,
    "minDist" = 0.1),
  generateNewClusterSettings = list(
    "reducedDims" = "Harmony",
    "method" = "Seurat",
    "name" = "Clusters2",
    "resolution" = 0.6
  ),
  filterProp = NA,
  runPreAnnot = FALSE,
  force = FALSE,
  runPreAnnotIsoComps = isoComparisonsART,
  prefixForGraphs = "art"
)

projART_matched <- artPreProcessed$archrProj
adtART_matched <- artPreProcessed$adtSeu
ARTAdtToUse <- artPreProcessed$adtToUse

rm(artPreProcessed)
gc()

save.image("allMatchedProjs.RData")

###############################################################################
# iv manual annot
###############################################################################
projInvitro_matched <- assignManualAnnotation(archrProj = projInvitro_matched,
  cluster = "Clusters2",
  annotFn = "manualClusterAnnotations/invitro.csv")

adtInvitro_matched$manualClusterAnnot <- projInvitro_matched$manualClusterAnnot
adtInvitro_matched$haystackOut <- projInvitro_matched$haystackOut

plotUmap(projInvitro_matched,
  colorBy = "manualClusterAnnot",
  fn = "InVitro_umap_labeledCluster",
  embedding = "UMAP2")
plotUmap(projInvitro_matched,
  colorBy = "Clusters2",
  colorLabelBy = "Clusters2",
  fn = "InVitro_umap_unannotCluster_withPlotLabels",
  embedding = "UMAP2",
  colorScheme = scale_color_manual(values = multiHueColorPalette),
  propInLegend = TRUE)
plotUmap(projInvitro_matched,
  colorBy = c("Clusters2", "manualClusterAnnot"),
  colorLabelBy = "Clusters2",
  fn = "InVitro_umap_labeledCluster_withPlotLabels",
  embedding = "UMAP2",
  customClusterSort = TRUE,
  colorScheme = scale_color_manual(values = multiHueColorPalette),
  propInLegend = TRUE)

plotDiscreteLollipop(projInvitro_matched, "inVitro_discreteAbsolute_matched",
  cluster = "manualClusterAnnot",
  graphType = "absolute")

plotDiscreteLollipop(projInvitro_matched, "inVitro_discreteHivOnly_matched",
  cluster = "manualClusterAnnot",
  graphType = "hivOnly")

###############################################################################
# art manual annot
###############################################################################

tmp <- loadArchRProject("../asapseq-hiv-art/A08A01B45A09_qcfiltTSS8_matched/")

tmpDf <- data.frame(
  cbc = projART_matched$cellNames,
  newAnnot = projART_matched$Clusters2)

tmpDf2 <- data.frame(
  cbc = tmp$cellNames,
  oldAnnot = tmp$Clusters2) %>%
  mutate(cbc = str_replace(cbc, "((A01|A08))_(\\d)", "\\1_post_\\2"))

tmpDf <- tmpDf %>%
  left_join(tmpDf2, by = "cbc") %>%
  dplyr::group_by(newAnnot, oldAnnot) %>%
  dplyr::tally() %>%
  group_by(newAnnot) %>%
  mutate(groupSum = sum(n)) %>%
  mutate(freqOfNew = n / groupSum)

ggplot(tmpDf, aes(x = newAnnot, y = oldAnnot, color = freqOfNew)) +
  geom_point() +
  theme_bw() +
  scale_color_distiller(palette = "YlOrRd", direction = 1)

t(adtART_matched@assays$tsa@data[c("A0034", "A0072", "A0154", "A0386", "A0390"), ]) %>%
  as.data.frame(.) %>%
  mutate(x = Cells(adtART_matched),
    y = adtART_matched$atacClusters) %>%
  tidyr::pivot_longer(cols = -one_of(c("x", "y")), names_to = "metric", values_to = "val") %>%
  group_by(metric, y) %>%
  dplyr::summarize(mean = mean(val)) %>%
  dplyr::arrange(metric, desc(y)) %>%
  tidyr::pivot_wider(names_from = "metric", values_from = "mean") %>%
  View(.)


projART_matched <- assignManualAnnotation(archrProj = projART_matched,
  cluster = "Clusters2",
  annotFn = "manualClusterAnnotations/art.csv")

adtART_matched$manualClusterAnnot <- projART_matched$manualClusterAnnot
adtART_matched$haystackOut <- projART_matched$haystackOut

plotUmap(projART_matched,
  colorBy = "manualClusterAnnot",
  fn = "art_umap_labeledCluster",
  embedding = "UMAP2")
plotUmap(projART_matched,
  colorBy = "Clusters2",
  colorLabelBy = "Clusters2",
  fn = "art_umap_unannotCluster_withPlotLabels",
  embedding = "UMAP2",
  colorScheme = scale_color_manual(values = multiHueColorPalette),
  propInLegend = TRUE)
plotUmap(projART_matched,
  fn = "art_umap_labeledCluster_withPlotLabels",
  colorBy = c("Clusters2", "manualClusterAnnot"),
  colorLabelBy = NULL,
  embedding = "UMAP2",
  customClusterSort = TRUE,
  colorScheme = scale_color_manual(values = multiHueColorPalette),
  propInLegend = TRUE)

plotDiscreteLollipop(projART_matched, "art_discreteAbsolute_matched",
  cluster = "manualClusterAnnot",
  gheight = 4,
  graphType = "absolute")

plotDiscreteLollipop(projART_matched, "art_discreteHivOnly_matched",
  cluster = "manualClusterAnnot",
  gheight = 4,
  graphType = "hivOnly")


artSampleOrder <- c("A01_pre", "A01_post", "A08_pre", "A08_post", "A09_pre", "A09_post", "B45")
plotDiscreteLollipop(projART_matched, "art_discreteAbsolute_matched_byDonor",
  cluster = "manualClusterAnnot",
  donorColumn = "individual",
  donorColumnOrder = artSampleOrder,
  gheight = 4.25,
  gwidth = 18,
  graphType = "absolute")

plotDiscreteLollipop(projART_matched, "art_discreteHivOnly_matched_byDonor",
  cluster = "manualClusterAnnot",
  donorColumn = "individual",
  donorColumnOrder = artSampleOrder,
  gheight = 4.25,
  gwidth = 16,
  showAggregate = TRUE,
  graphType = "hivOnly")


tmp <- generateUmapDfFromArchR(projART_matched,
  cluster = "manualClusterAnnot",
  secondGraphColumn = "haystackOut",
  donorColumn = "individual",
  customSort = TRUE)

tmp <- tmp %>%
  mutate(hivPos = ifelse(secondMetadata, "Pos", "Neg")) %>%
  dplyr::filter(hivPos != "Neg") %>%
  mutate(donor = factor(donor, levels = artSampleOrder))

tmp <- tmp %>% group_by(donor) %>% dplyr::count(cluster)

tmpAggr <- tmp %>%
  ungroup() %>%
  group_by(cluster) %>%
  dplyr::summarize(n = sum(n)) %>%
  mutate(donor = "Aggregate")

tmp <- tmp %>% 
  bind_rows(tmpAggr) %>%
  mutate(donor = factor(donor, levels = c(levels(tmp$donor), "Aggregate")))

tmp <- tmp %>%
  separate(donor, into = c("individual", "stage"), sep = "_", remove = FALSE) %>%
  mutate(ifelse(is.na(stage), "N/A", stage)) %>%
  mutate(stage = factor(stage, levels = c("pre", "post", "N/A"), labels = c("Pre", "Post", "N/A"))) %>%
  mutate(individual = factor(individual, levels = c("A01", "A08", "A09", "B45", "Aggregate"))) %>%
  ungroup() %>%
  group_by(donor) %>%
  mutate(proportion = n / sum(n) * 100)

ggplot(tmp, aes(x = proportion, y = cluster, fill = stage)) +
  geom_point(alpha = 0.5, color = "#000000", shape = 21, size = 2.75) +
  scale_y_discrete(limits = rev, drop = FALSE) +
  facet_wrap(~ individual, nrow = 1) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    text = element_text(family = "Arial", size = BASEPTFONTSIZE, color = "#000000"),
    axis.text = element_text(family = "Arial", size = BASEPTFONTSIZE, color = "#000000"),
    panel.grid.major.x = element_blank(),
    legend.text = element_text(family = "Arial", size = BASEPTFONTSIZE, color = "#000000"),
    axis.title = element_text(family = "Arial", size = BASEPTFONTSIZE, color = "#000000"),
    strip.background = element_rect(fill = NA),
    strip.text = element_text(family = "Arial", size = BASEPTFONTSIZE, color = "#000000")
  ) +
  labs(
    x = "% of sample (stage-specific except for Aggregate)",
    y = "Cluster",
    fill = "Stage")

###############################################################################
# iv differential analysis
###############################################################################
invitro_tcells_check <- grepl("^CD4", adtInvitro_matched$manualClusterAnnot)
invitro_tcells <- adtInvitro_matched$haystackOut[invitro_tcells_check]
invitro_markers <- findHIVDifferentialMarkers(seu = adtInvitro_matched,
  identPos = names(invitro_tcells[invitro_tcells]),
  identNeg = names(invitro_tcells[!invitro_tcells]),
  featuresToUse = invitroAdtToUse,
  tsa_catalog = tsa_catalog)

makeDifferentialLollipop(invitro_markers, "invitro_lollipop")
exportTsv(invitro_markers)

invitro_activeLaterCells <- adtInvitro_matched$haystackOut[grepl("CD4 Activated", adtInvitro_matched$manualClusterAnnot) & invitro_tcells_check]
invitro_activeLate_markers <- findHIVDifferentialMarkers(
  seu = adtInvitro_matched,
  tsa_catalog = tsa_catalog,
  identPos = names(invitro_activeLaterCells[invitro_activeLaterCells]),
  identNeg = names(invitro_activeLaterCells[!invitro_activeLaterCells]),
  featuresToUse = invitroAdtToUse)

makeDifferentialLollipop(invitro_activeLate_markers, "invitro_lollipop_activeLate")
exportTsv(invitro_activeLate_markers)

invitro_activeLate_markers_forVln <- invitro_activeLate_markers %>%
  group_by(Status) %>%
  arrange(desc(abs(piScore)), .by_group = TRUE) %>%
  top_n(10, abs(piScore))

plotVlnEnhanced(adtInvitro_matched, cells = names(invitro_activeLaterCells),
  feats = invitro_activeLate_markers_forVln$gene,
  xVar = "haystackOut",
  showAggregate = FALSE,
  separator = invitro_activeLate_markers_forVln$Status,
  titles = invitro_activeLate_markers_forVln$cleanName,
  fn = "invitro_activated_vln_top10")

invitro_earlyCells <- adtInvitro_matched$haystackOut[grepl("(Tcm|early)", adtInvitro_matched$manualClusterAnnot) & invitro_tcells_check]
invitro_early_markers <- findHIVDifferentialMarkers(
  seu = adtInvitro_matched,
  tsa_catalog = tsa_catalog,
  identPos = names(invitro_earlyCells[invitro_earlyCells]),
  identNeg = names(invitro_earlyCells[!invitro_earlyCells]),
  featuresToUse = invitroAdtToUse)

makeDifferentialLollipop(invitro_early_markers, "invitro_lollipop_early")
exportTsv(invitro_early_markers)

invitroGcPeakFn <- paste0(getOutputDirectory(projInvitro_matched), "_gcPeakCalled")
if (!dir.exists(invitroGcPeakFn)) {
  invitro_renamingDF <- getCellColData(projInvitro_matched, select = c("haystackOut", "manualClusterAnnot")) %>%
    as.data.frame(.) %>%
    mutate(cellNames = projInvitro_matched$cellNames) %>%
    mutate(condensedCluster = case_when(
      grepl("Activated", manualClusterAnnot) ~ "Activated",
      grepl("Tcm", manualClusterAnnot) ~ "Tcm",
      grepl("Naive", manualClusterAnnot) ~ "Naive_earlyMemory",
      grepl("Treg", manualClusterAnnot) ~ "Treg",
      TRUE ~ "Other")) %>%
    mutate(condensedHIVCluster = paste0(condensedCluster, "_", haystackOut))
  
  projInvitro_matched <- addCellColData(ArchRProj = projInvitro_matched,
    data = invitro_renamingDF$condensedHIVCluster,
    cells = invitro_renamingDF$cellNames,
    name = "condensedHivCluster",
    force = TRUE)
  
  projInvitro_matched <- addGroupCoverages(ArchRProj = projInvitro_matched,
    groupBy = "condensedHivCluster",
    useLabels = FALSE,
    minCells = 300,
    minReplicates = 3,
    force = TRUE)
  
  projInvitro_matched <- addReproduciblePeakSet(
    ArchRProj = projInvitro_matched,
    groupBy = "condensedHivCluster",
    peaksPerCell = 750,
    minCells = 100,
    pathToMacs2 = "/home/wuv/anaconda3/envs/asapseq/bin/macs2",
    force = TRUE
  )
  
  projInvitro_matched <- addPeakMatrix(projInvitro_matched, force = TRUE)
  projInvitro_matched <- addMotifAnnotations(ArchRProj = projInvitro_matched,
    motifSet = "cisbp",
    name = "Motif",
    force = TRUE)
  
  projInvitro_matched <- addBgdPeaks(projInvitro_matched)
  projInvitro_matched <- addDeviationsMatrix(
    ArchRProj = projInvitro_matched, 
    peakAnnotation = "Motif",
    force = TRUE
  )
  
  projInvitro_matched_gcPeak <- saveArchRProject(projInvitro_matched,
    outputDirectory = invitroGcPeakFn)
} else {
  projInvitro_matched_gcPeak <- loadArchRProject(invitroGcPeakFn, force = TRUE, showLogo = FALSE)
}

invitro_markerTest_Activated <- getMarkerFeatures(
  ArchRProj = projInvitro_matched_gcPeak, 
  useMatrix = "PeakMatrix",
  groupBy = "condensedHivCluster",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Activated_TRUE",
  bgdGroups = "Activated_FALSE",
  maxCells = 1000,
  verbose = FALSE
)

invitro_markerTest_activatedPeaks <- getMarkers(invitro_markerTest_Activated, cutOff = "FDR < 0.05", returnGR = TRUE)
invitro_markerTest_activatedPeaks <- invitro_markerTest_activatedPeaks$Activated_TRUE
invitro_markerTest_activatedPeaks <- hiAnnotator::getNearestFeature(invitro_markerTest_activatedPeaks,
  promoters(getGeneAnnotation(projInvitro_matched_gcPeak)[[1]], upstream = 0, downstream = 1), colnam = "nearestTSS",
  feature.colnam = "symbol")

invitro_markerTest_activatedPeaks <- hiAnnotator::getSitesInFeature(invitro_markerTest_activatedPeaks,
  getGeneAnnotation(projInvitro_matched_gcPeak)[[1]], colnam = "inGene",
  feature.colnam = "symbol")

invitro_markerTest_activatedPeaks_topPeaks <- data.frame(invitro_markerTest_activatedPeaks) %>%
  mutate(direction = ifelse(Log2FC > 0, "HIV+", "HIV-")) %>%
  mutate(piScore = Log2FC * -log10(FDR)) %>%
  group_by(direction) %>%
  arrange(abs(piScore), .by_group = TRUE) %>%
  ungroup() %>%
  mutate(symbol2 = case_when(
    inGene != "FALSE" & str_detect(as.character(inGene), nearestTSS) ~ glue("{nearestTSS}* ({nearestTSSDist})"),
    inGene != "FALSE" ~ glue("{inGene}*; {nearestTSS} ({nearestTSSDist})"),
    TRUE ~ glue("{nearestTSS} ({nearestTSSDist})"))) %>%
  mutate(symbol2 = str_wrap(symbol2, width = 35, exdent = 2)) %>%
  mutate(symbol2 = factor(symbol2, levels = symbol2))

exportTsv(invitro_markerTest_activatedPeaks_topPeaks)

iv_topPeaks_g <- ggplot(invitro_markerTest_activatedPeaks_topPeaks %>%
    filter(FDR < 0.05) %>%
    group_by(direction) %>%
    slice_max(order_by = abs(piScore), n = 15),
  aes(x = -log10(FDR), y = symbol2)) +
  geom_segment(aes(yend = symbol2, x = 0, xend = -log10(FDR)), linetype = "dashed", color = "#CCCCCC") +
  geom_point(aes(fill = Log2FC), color = "#000000", pch = 21, size = 2) +
  scale_fill_distiller(palette = "RdBu", direction = -1) +
  theme_classic() +
  theme(panel.background = element_rect(fill = NULL, colour = NULL),
    plot.margin = margin(0,0,0,0),
    strip.background = element_blank(),
    strip.text = element_text(size = BASEPTFONTSIZE, margin = margin(0, 0, 0, 0)),
    text = element_text(family = "Arial", size = BASEPTFONTSIZE),
    axis.text = element_text(size = 7, color = "#000000", family = "Arial"),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = BASEPTFONTSIZE),
    legend.box.margin = margin(0,0,0,-10),
    legend.direction = "vertical",
    legend.position = "right",
    legend.key.height = unit(0.5, "line"),
    legend.key.width = unit(0.5, "line")) +
  labs(fill = "log2(FC)") +
  scale_x_continuous(limits = c(0,NA)) +
  coord_cartesian(clip = "off") +
  scale_size_area(max_size = 4) +
  facet_wrap(~ direction, scales = "free")

savePlot(iv_topPeaks_g, fn = "invitro_markerTest_activated_topPeaks",
  devices = c("png", "rds"),
  gheight = 3.5, gwidth = 5)

plotVolcanoFromGetMarkerFeatures(invitro_markerTest_Activated,
  fn = "invitro_peaks_volcano_activatedHIV")

motifsUp <- peakAnnoEnrichment(
  seMarker = invitro_markerTest_Activated,
  ArchRProj = projInvitro_matched_gcPeak,
  peakAnnotation = "Motif",
  cutOff = "FDR < 0.05 & Log2FC > 0")

motifsDown <- peakAnnoEnrichment(
  seMarker = invitro_markerTest_Activated,
  ArchRProj = projInvitro_matched_gcPeak,
  peakAnnotation = "Motif",
  cutOff = "FDR < 0.05 & Log2FC < 0")

motifsAll <- peakAnnoEnrichment(
  seMarker = invitro_markerTest_Activated,
  ArchRProj = projInvitro_matched_gcPeak,
  peakAnnotation = "Motif",
  cutOff = "FDR < 0.05")

plotMotifRank(motifsDown, "invitro_motifsUp_activatedHIVneg", color = HIVNEGCOLOR)
plotMotifRank(motifsUp, "invitro_motifsUp_activatedHIVpos", color = HIVPOSCOLOR)

exportTsv(motifsAll, fn = "iv_activated_motifsAll", motifMode = TRUE)

iv_ccr5_track <- plotBrowserTrack2(projInvitro_matched_gcPeak,
  geneSymbol = "CCR5",
  useGroups = c("Activated_FALSE", "Activated_TRUE"),
  groupBy = "condensedHivCluster",
  upstream = 10000, downstream = 10000,
  scCellsMax = 750,
  plotSummary = c("bulkTrack", "scTracK", "geneTrack"),
  pal = c(HIVNEGCOLOR, HIVPOSCOLOR), borderWidth = 0, sizes = c(3,2,1))$CCR5

cleanUpTrackAndSave(iv_ccr5_track, fn = "invitro_activated_ccr5")

iv_genesOfInterest <- c("SELL", "IL7R", "PDCD1", "CD2", "ADGRG1")
for (tmp in iv_genesOfInterest) {
  cleanUpTrackAndSave(
    archrTrack = plotBrowserTrack2(projInvitro_matched_gcPeak,
      geneSymbol = tmp,
      useGroups = c("Activated_FALSE", "Activated_TRUE"),
      groupBy = "condensedHivCluster",
      upstream = 10000, downstream = 10000,
      scCellsMax = 750,
      plotSummary = c("bulkTrack", "scTracK", "geneTrack"),
      pal = c(HIVNEGCOLOR, HIVPOSCOLOR), borderWidth = 0, sizes = c(3,2,1))[[tmp]],
    fn = paste0("invitro_activated_", tmp))
} 

###############################################################################
# iv model analysis
###############################################################################
getAssayData <- function(seu, cellNames, adtToUse, splitRatio = 0.7, seed = 21) {
  set.seed(seed)
  
  df <- as.data.frame(t(seu@assays$tsa@data[adtToUse, cellNames])) %>%
    mutate(manualClusterAnnot = seu$manualClusterAnnot[cellNames],
      haystackOut = ifelse(seu$haystackOut[cellNames], 1, 0)) %>%
    dplyr::select(-manualClusterAnnot) %>%
    mutate(haystackOut = factor(haystackOut))
  
  sampleSplit <- caTools::sample.split(Y=df$haystackOut, SplitRatio = splitRatio) 
  trainSet <- subset(x = df, sampleSplit == TRUE)
  testSet <- subset(x = df, sampleSplit == FALSE)
  
  return(list("train" = trainSet, "test" = testSet))
}

calcLogisticalRegression <- function(trainSet) {
  lmModel <- glm(haystackOut ~ ., family = binomial(link='logit'), data = trainSet)
  
  return(lmModel)
}

calcNaiveBayes <- function(trainSet) {
  bayesModel <- naivebayes::naive_bayes(haystackOut ~ ., data = trainSet, useKernel = TRUE)
  return(bayesModel)
}

calcRandomForest <- function(trainSet, sampleProps = list(c(NULL, NULL), c(2,1), c(1, 1), c(0.5, 1))) {
  trainSetTable <- table(trainSet$haystackOut)
  
  res <- lapply(sampleProps, function(prop) {
    if (is.null(prop)) {
      propCalculated <- sum(trainSetTable)
    } else {
      propCalculated <- c(trainSetTable[[2]] * prop[[1]], trainSetTable[[2]] * prop[[2]])
    }
    
    rf <- randomForest::randomForest(
      haystackOut ~ ., trainSet,
      ntree = 4000,
      sampsize = propCalculated,
      strata = trainSet$haystackOut,
      importance = TRUE)
  })
  
  names(res) <- paste("RF (", sapply(sampleProps, paste, collapse = ":"), ")", sep = "")
  names(res)[grep("RF \\(\\)", names(res))] <- "RF (all)"
  
  return(res)
}

getClassifierPrediction <- function(predicted, actual) {
  if (class(predicted) == "prediction") {
    pred <- predicted
  } else {
    pred <- ROCR::prediction(predicted, actual)
  }
  
  return(pred)
}

getClassifierPerformance <- function(pred, perfMode = "roc") {
  if (perfMode == "roc") {
    perf <- ROCR::performance(pred, measure='tpr', x.measure='fpr')
    
  } else if (perfMode == "auc") {
    perf <- ROCR::performance(pred, measure = perfMode)
    
    return(perf@y.values[[1]])
    
  } else {
    perf <- ROCR::performance(pred, measure = perfMode)
  }
  
  return(perf)
}


invitroModelData <- getAssayData(adtInvitro_matched, names(invitro_tcells), invitroAdtToUse)
invitroLR <- calcLogisticalRegression(invitroModelData[["train"]])
invitroNB <- calcNaiveBayes(invitroModelData[["train"]])
invitroRF <- calcRandomForest(invitroModelData[["train"]])

invitroRFPred <- lapply(invitroRF, function(x) {
  pred <- getClassifierPrediction(predict(x, invitroModelData[["test"]], type = "prob")[, 2],
    invitroModelData[["test"]]$haystackOut)
  
  return(pred)
})
names(invitroRFPred) <- names(invitroRF)

invitroModelsPred <- list(
  "logistic" = getClassifierPrediction(predict(invitroLR, invitroModelData[["test"]], type = 'response'),
    invitroModelData[["test"]]$haystackOut),
  "naive bayes" = getClassifierPrediction(predict(invitroNB, invitroModelData[["test"]], type = "prob")[, 2],
    invitroModelData[["test"]]$haystackOut)
)

invitroModelsPred <- base::append(invitroModelsPred, invitroRFPred)
invitroModelsROC <- lapply(invitroModelsPred, getClassifierPerformance)
invitroModelsAUC <- lapply(invitroModelsPred, getClassifierPerformance, perfMode = "auc")

rm(invitroModelData)

plotGgRoc(invitroModelsROC, invitroModelsAUC, "invitro_models")

plotLogitRegressionVolcano(invitroLR, tsa_catalog, invitroAdtToUse, "invitro_logitRegression_volcano")

randomForest::importance(invitroRF$`RF (all)`, type = "1") %>% 
  as.data.frame(.) %>% 
  arrange(desc(MeanDecreaseAccuracy)) %>%
  mutate(DNA_ID = rownames(.)) %>%
  left_join(tsa_catalog %>% select(DNA_ID, cleanName), by = "DNA_ID") %>%
  head(n = 10)

# invitro activated
invitroActModelData <- getAssayData(adtInvitro_matched, names(invitro_activeLaterCells), invitroAdtToUse)
invitroActLR <- calcLogisticalRegression(invitroActModelData[["train"]])
invitroActNB <- calcNaiveBayes(invitroActModelData[["train"]])
invitroActRF <- calcRandomForest(invitroActModelData[["train"]], sampleProps = list(c(NULL, NULL), c(1, 1), c(0.5, 1)))

invitroActRFPred <- lapply(invitroActRF, function(x) {
  pred <- getClassifierPrediction(predict(x, invitroActModelData[["test"]], type = "prob")[, 2],
    invitroActModelData[["test"]]$haystackOut)
  
  return(pred)
})
names(invitroActRFPred) <- names(invitroActRF)

invitroActModelsPred <- list(
  "logistic" = getClassifierPrediction(predict(invitroActLR, invitroActModelData[["test"]], type = 'response'),
    invitroActModelData[["test"]]$haystackOut),
  "naive bayes" = getClassifierPrediction(predict(invitroActNB, invitroActModelData[["test"]], type = "prob")[, 2],
    invitroActModelData[["test"]]$haystackOut)
)

invitroActModelsPred <- base::append(invitroActModelsPred, invitroActRFPred)
invitroActModelsROC <- lapply(invitroActModelsPred, getClassifierPerformance)
invitroActModelsAUC <- lapply(invitroActModelsPred, getClassifierPerformance, perfMode = "auc")

rm(invitroActModelData)

plotGgRoc(invitroActModelsROC, invitroActModelsAUC,"invitro_models_activated")

randomForest::importance(invitroActRF$`RF (all)`, type = "1") %>% 
  as.data.frame(.) %>% 
  arrange(desc(MeanDecreaseAccuracy)) %>%
  mutate(DNA_ID = rownames(.)) %>%
  left_join(tsa_catalog %>% select(DNA_ID, cleanName), by = "DNA_ID") %>%
  head(n = 10)


###############################################################################
# art differential analysis
###############################################################################
art_tcells_check <- (grepl("CD4", adtART_matched$manualClusterAnnot))
art_tcells <- adtART_matched$haystackOut[art_tcells_check]
art_markersLR <- findHIVDifferentialMarkers(
  seu = adtART_matched,
  identPos = names(art_tcells[art_tcells]),
  identNeg = names(art_tcells[!art_tcells]),
  featuresToUse = ARTAdtToUse,
  findMarkerMethod = "LR",
  logfc.threshold = 0.1,
  tsa_catalog = tsa_catalog,
  latent.vars = c("individual", "manualClusterAnnot"))

makeDifferentialLollipop(art_markersLR, "art_lollipop_LR")
exportTsv(art_markersLR)

art_markersDEseq2 <- findHIVDifferentialMarkers(
  seu = adtART_matched,
  identPos = names(art_tcells[art_tcells]),
  identNeg = names(art_tcells[!art_tcells]),
  featuresToUse = ARTAdtToUse,
  findMarkerMethod = "DESeq2",
  logfc.threshold = 0.1,
  tsa_catalog = tsa_catalog)

makeDifferentialLollipop(art_markersDEseq2, "art_lollipop_deseq2")
exportTsv(art_markersDEseq2)

art_markersWilcox <- findHIVDifferentialMarkers(
  seu = adtART_matched,
  identPos = names(art_tcells[art_tcells]),
  identNeg = names(art_tcells[!art_tcells]),
  featuresToUse = ARTAdtToUse,
  logfc.threshold = 0.1,
  findMarkerMethod = "wilcox",
  tsa_catalog = tsa_catalog)

makeDifferentialLollipop(art_markersWilcox, "art_lollipop_wilcox")
exportTsv(art_markersWilcox)


plotVlnEnhanced(adtART_matched, cells = names(art_tcells),
  feats = art_markersWilcox$gene,
  separator = art_markersWilcox$Status,
  titles = art_markersWilcox$cleanName,
  fn = "art_tcells_vln_wilcox")

### run on a less condensed HIV Cluster grouping
artTestingDf <- data.frame(
  manual = projART_matched$manualClusterAnnot,
  cellNames = projART_matched$cellNames,
  haystackOut = projART_matched$haystackOut) %>%
  mutate(lessCondensed = case_when(
    !grepl("^CD4", manual) ~ "Other",
    grepl("MAIT", manual) ~ "MAIT",
    grepl("(Tcm|cTfh)", manual) ~ "Tcm/Ttm",
    TRUE ~ "Tem/effector"
  )) %>%
  mutate(lessCondensedHivClust = paste(lessCondensed, haystackOut, sep = "_"))

projART_matched <- addCellColData(ArchRProj = projART_matched,
  data = artTestingDf$lessCondensedHivClust,
  cells = artTestingDf$cellNames,
  name = "lessCondensedHivClust",
  force = TRUE)

rm(artTestingDf)

artLessGcPeakFn <- paste0(getOutputDirectory(projART_matched), "_lessCondensed_gcPeakCalled")
if (!dir.exists(artLessGcPeakFn)) {
  projART_matchedLess_subset <- projART_matched[!grepl("(Other)", projART_matched$lessCondensedHivClust), ]
  
  projART_matchedLess_subset <- addGroupCoverages(
    ArchRProj = projART_matchedLess_subset,
    groupBy = "lessCondensedHivClust",
    minReplicates = 2,
    minCells = 29,
    force = TRUE,
    sampleRatio = 0.9)
  
  projART_matchedLess_subset <- addReproduciblePeakSet(
    ArchRProj = projART_matchedLess_subset,
    groupBy = "lessCondensedHivClust", 
    pathToMacs2 = "/home/wuv/anaconda3/envs/asapseq/bin/macs2")
  
  projART_matchedLess_subset <- addPeakMatrix(projART_matchedLess_subset, force = TRUE)
  projART_matchedLess_subset <- addMotifAnnotations(ArchRProj = projART_matchedLess_subset, motifSet = "cisbp", name = "Motif", force = TRUE)
  projART_matchedLess_subset <- addBgdPeaks(projART_matchedLess_subset, force = TRUE)
  
  projART_matchedLess_subset <- addDeviationsMatrix(
    ArchRProj = projART_matchedLess_subset, 
    peakAnnotation = "Motif",
    force = TRUE
  )
  
  projART_matchedLess_subset_gcPeak <- saveArchRProject(
    projART_matchedLess_subset,
    outputDirectory = artLessGcPeakFn,
    load = TRUE,
    dropCells = TRUE)
  
} else {
  projART_matchedLess_subset_gcPeak <- loadArchRProject(artLessGcPeakFn)
}


##############
### Checkpoint
##############
save.image(paste(format(Sys.time(), "%Y-%m-%d-%H-%M-rebuttal"), "RData", sep = "."))


artLess_markerTest_deviations <- getMarkerFeatures(
  ArchRProj = projART_matchedLess_subset_gcPeak, 
  useMatrix = "MotifMatrix",
  groupBy = "lessCondensedHivClust",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  verbose = FALSE,
  useSeqnames = "z",
)

# plots the mean z score for each motif...

heatmapCutoff <- "FDR < 0.05 & MeanDiff > 0.5"
artLess_markerTest_deviations_heatmap <- plotMarkerHeatmap(
  seMarker = artLess_markerTest_deviations, 
  log2Norm = FALSE,
  cutOff = heatmapCutoff,
  transpose = FALSE,
  returnMatrix = TRUE,
  scaleRows = FALSE
)

artLess_markerTest_deviations_heatmapText <- bind_rows(
  lapply(getMarkers(artLess_markerTest_deviations, cutOff = heatmapCutoff), function(x) {
    return(as.data.frame(x))
  }
  ), .id = "group")

artLess_markerTest_deviations_heatmapSig <- matrix(" ",
  nrow = nrow(artLess_markerTest_deviations_heatmap),
  ncol = ncol(artLess_markerTest_deviations_heatmap))

dimnames(artLess_markerTest_deviations_heatmapSig) <- dimnames(artLess_markerTest_deviations_heatmap)

for (i in seq_along(artLess_markerTest_deviations_heatmapText$group)) {
  row <- artLess_markerTest_deviations_heatmapText[i, ]
  
  artLess_markerTest_deviations_heatmapSig[row$name, row$group] <- case_when(
    row$FDR < 0.01 ~ "✱✱",
    row$FDR < 0.05 ~ "✱",
    TRUE ~ " ")
}

saveRDS(artLess_markerTest_deviations_heatmap, file = "outs/rds/artLess_markerTest_deviations_heatmap.rds")
saveRDS(artLess_markerTest_deviations_heatmapSig, file = "outs/rds/artLess_markerTest_deviations_heatmapTxt.rds")
tmp <- ComplexHeatmap::Heatmap(
  matrix = artLess_markerTest_deviations_heatmap,
  name = " ",
  border = FALSE,
  cluster_rows = TRUE,
  show_row_names = FALSE,
  row_dend_side = "left",
  show_row_dend = TRUE,
  cluster_columns = FALSE,
  column_split = c(1,1,2,2,3,3),
  show_column_names = FALSE,
  show_heatmap_legend = TRUE,
  column_dend_side = "top",
  column_title = NULL,
  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
    grid.text(artLess_markerTest_deviations_heatmapSig[i, j], x, y, gp = gpar(fontsize = 2.5, hjust = 0.5, vjust = 0.5))
  },
  heatmap_legend_param = list(
    legend_direction = "horizontal",
    legend_width = unit(3, "cm")),
  show_column_dend = TRUE,
  # column_km = 3,
  row_km = 10,
  bottom_annotation = ComplexHeatmap::HeatmapAnnotation(
    `Phenotype` = str_split_fixed(colnames(artLess_markerTest_deviations_heatmap), "_", n = 2)[, 1],
    `HIV` = ifelse(str_split_fixed(colnames(artLess_markerTest_deviations_heatmap), "_", n = 2)[, 2] == "TRUE", "HIV+", "HIV-"),
    col = list(
      `Phenotype` = c("Tem/effector" = "#1b9e77", "MAIT" = "#d95f02", "Tcm/Ttm" = "#7570b3"),
      `HIV` = c("HIV+" = HIVPOSCOLOR, "HIV-" = HIVNEGCOLOR)
    ),
    annotation_name_gp = gpar(fontsize = 6)
  ),
  right_annotation = rowAnnotation(
    foo = anno_block(gp = gpar(fill = "#000000"), width = unit(1, "mm")),
    bar = anno_block(
      graphics = function(index, levels) {
        lbls <- rownames(artLess_markerTest_deviations_heatmap)[index]
        lbls <- gsub("_\\d+", "", lbls)
        lbls <- str_wrap(paste(lbls, collapse = " "), width = 40)
        
        grid.rect(gp = gpar(fill = NA, col = NA))
        txt = paste(lbls, collapse = ",")
        grid.text(txt, 0.01, 0.5, rot = 0, hjust = 0, gp = gpar(fontsize = 8))
      },
      width = unit(7, "cm"))
  )
)


artLess_markerTest_tem_deviations <- getMarkerFeatures(
  ArchRProj = projART_matchedLess_subset_gcPeak, 
  useMatrix = "MotifMatrix",
  useGroups = "Tem/effector_TRUE",
  bgdGroups = "Tem/effector_FALSE",
  closest = FALSE,
  groupBy = "lessCondensedHivClust",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  verbose = FALSE,
  useSeqnames = "z",
)
plotMotifDot(artLess_markerTest_tem_deviations, "art_chromVAR_temEffector_motifsUp_HIVneg", direction = "negative", pValMax = 0.05, showTopNByPVal = 10)
plotMotifDot(artLess_markerTest_tem_deviations, "art_chromVAR_temEffector_motifsUp_HIVpos", direction = "positive", pValMax = 0.05, showTopNByPVal = 10)

artLess_markerTest_tcmttm_deviations <- getMarkerFeatures(
  ArchRProj = projART_matchedLess_subset_gcPeak, 
  useMatrix = "MotifMatrix",
  useGroups = "Tcm/Ttm_TRUE",
  bgdGroups = "Tcm/Ttm_FALSE",
  closest = FALSE,
  groupBy = "lessCondensedHivClust",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  verbose = FALSE,
  useSeqnames = "z",
)
# plotMotifDot(artLess_markerTest_tcmttm_deviations,
#   "art_chromVAR_tcmttm_motifsUp_HIVneg", direction = "negative", pValMax = 0.05, showTopNByPVal = 15, alphaBySig = TRUE)
# plotMotifDot(artLess_markerTest_tcmttm_deviations,
#   "art_chromVAR_tcmttm_motifsUp_HIVpos", direction = "positive", pValMax = 0.05, showTopNByPVal = 15, alphaBySig = TRUE)

plotEnhancedMotifDot(
  markerFeatures = artLess_markerTest_tcmttm_deviations,
  fn = "art_chromVAR_tcmttm_motifsUp_HIVpos_enhanced",
  devices = "svg",
  preFilter = "meandiff > 0",
  xVar = "motif",
  yVar = "fdr",
  subtitle = "Tcm/Ttm HIV+ enriched",
  colorVar = "meandiff",
  showTopNByPVal = 20
)

plotEnhancedMotifDot(
  markerFeatures = artLess_markerTest_tcmttm_deviations,
  fn = "art_chromVAR_tcmttm_motifsUp_HIVneg_enhanced",
  preFilter = "meandiff < 0",
  xVar = "motif",
  yVar = "fdr",
  subtitle = "Tcm/Ttm HIV- enriched",
  colorVar = "meandiff",
  showTopNByPVal = 15
)

# artLess_markerTest_mait_deviations <- getMarkerFeatures(
#   ArchRProj = projART_matchedLess_subset_gcPeak, 
#   useMatrix = "MotifMatrix",
#   useGroups = "MAIT_TRUE",
#   bgdGroups = "MAIT_FALSE",
#   closest = FALSE,
#   groupBy = "lessCondensedHivClust",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   verbose = FALSE,
#   useSeqnames = "z",
# )
# plotMotifDot(artLess_markerTest_mait_deviations,
#   "art_chromVAR_mait_motifsUp_HIVneg", direction = "negative", pValMax = 0.05, showTopNByPVal = 10)
# plotMotifDot(artLess_markerTest_mait_deviations,
#   "art_chromVAR_mait_motifsUp_HIVpos", direction = "positive", pValMax = 0.05, showTopNByPVal = 10)

# adt analysis on the less condensed cluster form
adtART_matched$lessCondensedHivClust <- projART_matched$lessCondensedHivClust

art_tcmttm_markers <- findHIVDifferentialMarkers(seu = adtART_matched,
  identPos = names(adtART_matched$lessCondensedHivClust[adtART_matched$lessCondensedHivClust == "Tcm/Ttm_TRUE" & adtART_matched$manualClusterAnnot != "CD4 Tcm/Ttm (recently activated)"]),
  identNeg = names(adtART_matched$lessCondensedHivClust[adtART_matched$lessCondensedHivClust == "Tcm/Ttm_FALSE" & adtART_matched$manualClusterAnnot != "CD4 Tcm/Ttm (recently activated)"]),
  featuresToUse = ARTAdtToUse,
  findMarkerMethod = "DESeq2",
  logfc.threshold = 0.1,
  tsa_catalog = tsa_catalog)

makeDifferentialLollipop(art_tcmttm_markers, "art_tcmttm_markers")
exportTsv(art_tcmttm_markers)

art_tem_markers <- findHIVDifferentialMarkers(seu = adtART_matched,
  identPos = names(adtART_matched$lessCondensedHivClust[adtART_matched$lessCondensedHivClust == "Tem/effector_TRUE"]),
  identNeg = names(adtART_matched$lessCondensedHivClust[adtART_matched$lessCondensedHivClust == "Tem/effector_FALSE"]),
  featuresToUse = ARTAdtToUse,
  findMarkerMethod = "DESeq2",
  logfc.threshold = 0.1,
  tsa_catalog = tsa_catalog)

makeDifferentialLollipop(art_tem_markers, "art_tem_markers")
exportTsv(art_tem_markers)

# art_mait_markers <- findHIVDifferentialMarkers(seu = adtART_matched,
#   identPos = names(adtART_matched$lessCondensedHivClust[adtART_matched$lessCondensedHivClust == "MAIT_TRUE"]),
#   identNeg = names(adtART_matched$lessCondensedHivClust[adtART_matched$lessCondensedHivClust == "MAIT_FALSE"]),
#   featuresToUse = ARTAdtToUse,
#   findMarkerMethod = "wilcox",
#   logfc.threshold = 0.1,
#   tsa_catalog = tsa_catalog)
# 
# makeDifferentialLollipop(art_mait_markers, "art_mait_markers")
# exportTsv(art_mait_markers)

###############################################################################
# qc and frags
###############################################################################
plotUpsetCBC(atacCBC = projInVitro$cellNames,
  hivCBC = unique(haystackInVitro$allViralFrags$newCbc),
  adtCBC = Cells(adtInVitro),
  fn = "invitro_upset_cbc")

plotUpsetCBC(atacCBC = projART$cellNames,
  hivCBC = unique(haystackART$allViralFrags$newCbc),
  adtCBC = Cells(adtART),
  separateByIndividual = TRUE,
  gheight = 4,
  fn = "art_upset_cbc")

calculateLTRVersusInternal <- function(frags, ltr5End, ltr3Start, ltr3End) {
  df <- frags %>%
    mutate(ltrOnly = (startBp <= ltr5End & endBp <= ltr5End) | (startBp >= ltr3Start & endBp <= ltr3End)) %>%
    mutate(internalOnly = startBp > ltr5End & endBp < ltr3Start) %>%
    mutate(both = !ltrOnly & !internalOnly) %>%
    summarize(LTRonly = sum(ltrOnly),
      internalOnly = sum(internalOnly),
      both = sum(both)) %>%
    pivot_longer(everything(), names_to = "metric", values_to = "n") %>%
    mutate(prop = round(n / sum(n), digits = 2))
  
  return(df)
}

sumaGeneAnnots <- read.csv("viralGenomeAnnotations/suma.csv") %>%
  mutate(annotation = factor(annotation, levels = rev(unique(annotation)))) %>%
  mutate(annotation2 = ifelse(annotation == "LTR", paste0(annotation, row_number()), annotation))

inVitroFrags_matched <- haystackInVitro$filteredviralFrags[haystackInVitro$filteredviralFrags$newCbc %in% projInvitro_matched$cellNames, ]

plotFragMultiGraph(
  frags = inVitroFrags_matched, coverage = rep(0, 9726),
  geneAnnots = sumaGeneAnnots, fn = "invitro_frags")
calculateLTRVersusInternal(inVitroFrags_matched, 632, 9093, 9725)

hxb2GeneAnnots <- read.csv("viralGenomeAnnotations/hxb2.csv") %>%
  mutate(annotation = factor(annotation, levels = rev(unique(annotation)))) %>%
  mutate(annotation2 = ifelse(annotation == "LTR", paste0(annotation, row_number()), annotation))

artFrags_matched <- haystackART$filteredviralFrags[haystackART$filteredviralFrags$newCbc %in% projART_matched$cellNames, ]
artFrags_matched <- artFrags_matched %>%
  mutate(id = paste0(sample, cbc, seqname, startBp, endBp, readname, sep = "")) %>%
  mutate(individual = str_match(sample, "(.*)(?=_.*$)")[, 1])

checkFragmentOverlap <- function(q1, q2, s1, s2) {
  if (q1 <= s2 & q1 >= s1) {
    return(TRUE)
  } else if (q2 >= s1 & q2 <= s2) {
    return(TRUE)
  } else if (s2 >= q1 & s2 <= q2) {
    return(TRUE)
  } else if (s1 >= q1 & s1 <= q2) {
    return(TRUE)
  }
  return(FALSE)
}

artConsensusAnnots <- read.csv("viralGenomeAnnotations/GeneCutterParser_GeneCutterParser_A08_A01_B45_A09.tsv", sep = "\t") %>%
  filter(annotation != "Genome") %>%
  mutate(annotation = tolower(annotation)) %>%
  mutate(annotation = ifelse(grepl("ltr", annotation), "LTR", annotation)) %>%
  bind_rows(hxb2GeneAnnots %>% dplyr::select(-annotation2) %>% mutate(genome = "chrHXB2")) %>%
  group_split(genome)

names(artConsensusAnnots) <- sapply(artConsensusAnnots, function(x) {return(x[1, "genome"])})

artFrags_matchedWithAnnot <- lapply(seq_along(artFrags_matched$sample), function(i) {
  seqname <- artFrags_matched[i, "seqname"]
  seqnameAnnots <- artConsensusAnnots[[seqname]]
  
  if (is.null(seqnameAnnots)) {
    return(NULL)
  }
  
  q1 <- artFrags_matched[i, "startBp"]
  q2 <- artFrags_matched[i, "endBp"]
  
  id <- paste(artFrags_matched[i, "sample"],
    artFrags_matched[i, "cbc"],
    artFrags_matched[i, "seqname"],
    artFrags_matched[i, "startBp"],
    artFrags_matched[i, "endBp"],
    artFrags_matched[i, "readname"], sep = "")
  
  matches <- lapply(seq_along(seqnameAnnots$annotation), function(i) {
    s1 <- seqnameAnnots$startPos[i]
    s2 <- seqnameAnnots$endPos[i]
    
    check <- checkFragmentOverlap(q1, q2, s1, s2)
    
    if (check) {
      return(data.frame(id = id, annot = seqnameAnnots$annotation[i]))
    } else {
      return(NULL)
    }
  })
  
  return(bind_rows(matches))
})


artFrags_matchedWithAnnot <- bind_rows(artFrags_matchedWithAnnot)
plogFragMultiAnnotGraph(artFrags_matched, artFrags_matchedWithAnnot, unique(hxb2GeneAnnots$annotation), fn = "art_frags2")

uniqArtSamples <- unique(artFrags_matched$individual)
for (smpl in uniqArtSamples) {
  plogFragMultiAnnotGraph(
    artFrags_matched %>% filter(individual == smpl),
    artFrags_matchedWithAnnot,
    unique(hxb2GeneAnnots$annotation),
    fn = paste0("art_frags2_", smpl))
  
}

makeFancyUpsetPlotHIV(
  seu = adtInvitro_matched,
  cellsOfInterest = names(invitro_activeLaterCells),
  featuresOfInterest = c("A0141", "A0870", "A0367", "A0147", "A0390"),
  featuresOfInterestNames = c("CCR5", "SLAM", "CD2", "CD62L", "CD127"),
  metadata = ifelse(invitro_activeLaterCells, "HIV+", "HIV-"),
  thresholds = c(1.25, 1.25, 1.25, 2.75, 1.5),
  fn = "invitro_fancyUpset_activatedLaterCells"
)

save.image()

salantesData <- data.frame(
  sample = c("A09_pre", "A09_post", "A01", "A08"),
  iupm = c(0.614, 1.237, 8.35, 0.141),
  total_dna = c(1208.3, 1221.8, 293.8, 1564.5),
  rna = c(10834.7, 7203.2, 754.1, 4528.3)
)

rownames(salantesData) <- salantesData$sample

data.frame(table(projART_matched$individual, projART_matched$haystackOut)) %>%
  as_tibble() %>%
  dplyr::rename(sample = Var1,
    hivPos = Var2,
    freq = Freq) %>%
  group_by(sample) %>%
  mutate(sample_total = sum(freq)) %>%
  group_by(sample, hivPos) %>%
  mutate(prop = freq / sample_total * 100) %>%
  mutate(sample = as.character(sample)) %>%
  mutate(hivPos = as.logical(hivPos)) %>%
  filter(hivPos & sample != "B45") %>% 
  left_join(salantesData) %>%
  select(sample, hivPos, 5:8) %>%
  mutate(total_dna = log10(total_dna)) %>%
  mutate(rna = log10(rna)) %>%
  mutate(iupm = log10(iupm)) %>%
  pivot_longer(iupm:rna, names_to = "metric", values_to = "value") %>%
  ggplot(aes(y = prop, x = value)) +
  geom_point(aes(color = sample)) +
  facet_wrap(~ metric, scales = "free_x") +
  geom_smooth(method='lm', formula= y ~ x) +
  stat_cor(label.y = 0.3, method = "spearman", cor.coef.name = "rho") +
  theme_bw() +
  theme(panel.grid = element_blank())


tmp2 <- invitro_markerTest_activatedPeaks_topPeaks %>%
  mutate(finder = gsub(",", "|", inGene)) %>%
  mutate(finder = paste0("(", finder, ")")) %>%
  mutate(finder2 = gsub(",", "|", nearestTSS)) %>%
  mutate(finder2 = paste0("(", finder, ")"))


tmp2$inGeneInteractHIV <- sapply(tmp2$finder, function(x) {
  hit <- grepl(x, tmp$Human_GeneSymbol)
  return(paste(tmp$Interaction_Desc[hit], collapse = " | "))
  
})
  
tmp2$nearestInteractHIV <- sapply(tmp2$finder2, function(x) {
  hit <- grepl(x, tmp$Human_GeneSymbol)
  return(paste(tmp$Interaction_Desc[hit], collapse = " | "))
})


# makeFancyUpsetPlotHIV(
#   seu = adtART_matched,
#   cellsOfInterest = names(art_tcells),
#   featuresOfInterest = c("A0047", "A0576", "A0396", "A0575", "A0088"),
#   featuresOfInterestNames = c("CD200 (OX2)", "CD49d", "CD26", "CD49a", "PD-1"),
#   metadata = ifelse(art_tcells, "HIV+", "HIV-"),
#   thresholds = c(1, 1.75, 2, 1, 1, 1.5),
#   fn = "art_fancyUpset_tCells"
# )

# tmp <- calculateCoverage(inVitroFrags_matched, rep(0, 9726))
# 
# tmp2 <- calculateCoverage(inVitroFrags_matched %>%
#     group_by(readname, newCbc) %>%
#     filter(n() == 2) %>%
#     arrange(min(startBp, endBp), .by_group = TRUE) %>%
#     summarise(inferredStart = first(endBp) + 1, inferredEnd = last(startBp) -1) %>%
#     as.data.frame(.),
#   coverage = rep(0, 9726), startCol = "inferredStart", endCol = "inferredEnd")
# 
# tmp %>%
#   mutate(y = y + tmp2$y) %>%
#   filter(x < 800) %>%
#   ggplot(aes(x = x, y=y)) +
#   annotate("rect",
#     xmin = c(40, 456-2), xmax = c(200, 456+140), ymin = 0, ymax = 0.3,
#     alpha = .1, fill = "blue") +
#   annotate("rect",
#     xmin = c(315), xmax = c(408), ymin = 0, ymax = 0.3,
#     alpha = .1, fill = "red") +
#   geom_vline(xintercept = 456, color = "red") +
#   geom_line() +
#   theme_classic()


# mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# tmpSapEnsembl <- tsa_catalog$Ensembl.ID[tsa_catalog$Ensembl.ID != ""]
# names(tmpSapEnsembl) <- tsa_catalog[tsa_catalog$Ensembl.ID != "", "DNA_ID"]
# 
# results <- getBM(
#   attributes = c("ensembl_gene_id", "go_id","name_1006"),
#   filters = "ensembl_gene_id",
#   values = tmpSapEnsembl,
#   mart = mart)
# 
# ignoredGoTerms <- c(
#   "integral component of plasma membrane",
#   "membrane",
#   "plasma membrane",
#   "integral component of membrane",
#   "external side of plasma membrane",
#   ""
# )
# 
# tmpCleanPrefix <- function(x) {
#   return(gsub("(.+___)","",x))
# }
# 
# art_markers %>%
#   select(gene, cleanName, piScore, Status) %>%
#   mutate(geneID = tmpSapEnsembl[gene]) %>%
#   left_join(results, by = c("geneID" = "ensembl_gene_id")) %>%
#   group_by(Status, name_1006) %>%
#   tally() %>%
#   filter(!name_1006 %in% ignoredGoTerms) %>%
#   filter(n >= 2) %>%
#   arrange(n, .by_group = TRUE) %>%
#   mutate(xvar = paste0(Status, "___", name_1006)) %>%
#   mutate(xvar = factor(xvar, levels = unique(xvar))) %>%
#   slice_max(n, n = 20) %>%
#   ggplot(aes(x = n, y = xvar, fill = Status)) +
#   geom_col(color = "#000000") +
#   theme(
#     legend.position = "none",
#     panel.background = element_blank(), 
#     panel.grid.major.x = element_line(color = "#cccccc", linetype = "dotted"),
#     axis.text.y = element_text(size = 8),
#     axis.title = element_blank()) +
#   scale_y_discrete(labels = tmpCleanPrefix) +
#   scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
#   scale_fill_manual(values = c(HIVNEGCOLOR, HIVPOSCOLOR)) +
#   facet_grid(Status ~ ., scales = "free", space='free')
  

