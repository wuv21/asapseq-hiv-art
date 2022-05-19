###############################################################################
# pkg and custom scripts
###############################################################################
suppressMessages(library(tidyverse))
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
load("rds/chronic_seuratMerged.RData")
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
# load in chronic project from preHPC 
###############################################################################
projChronic <- loadArchRProject(path = "C01C02_qcfiltTSS6/", showLogo = FALSE)

uniqChronicSamples <- unique(projChronic$Sample)
haystackSamplePairingChronic <- paste0(uniqChronicSamples, "_hxb2")
names(haystackSamplePairingChronic) <- uniqChronicSamples

haystackChronic <- addHaystackData(proj = projChronic,
  haystackParentDir = "data/hiv-haystack",
  haystackSamples = haystackSamplePairingChronic)

projChronic <- haystackChronic$newProj

###############################################################################
# load in art project from preHPC 
###############################################################################
projART <- loadArchRProject(path = "A08A01B45A09_qcfiltTSS8/", showLogo = FALSE)

# add haystack data
uniqSamples <- unique(projART$Sample)
haystackSamplePairingART <- sapply(uniqSamples, function(i) {
  if (grepl("A01", i)) {
    newSampleName <- paste0(i, "_v3WithEnv_filtered")
  } else if (grepl("A08", i)) {
    newSampleName <- paste0(i, "_v3")
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
    "nNeighbors" = 200,
    "force" = TRUE,
    "minDist" = 0.3),
  generateNewClusterSettings = list(
    "method" = "Seurat",
    "name" = "Clusters2",
    "resolution" = 0.8
  ),
  runPreAnnot = TRUE,
  force = TRUE,
  runPreAnnotIsoComps = isoComparisonsInvitro,
  prefixForGraphs = "invitro"
)

projInvitro_matched <- ivPreProcessed$archrProj
adtInvitro_matched <- ivPreProcessed$adtSeu
invitroAdtToUse <- ivPreProcessed$adtToUse

rm(ivPreProcessed)

###############################################################################
# chronic processing
###############################################################################

chronicPreProcessed <- loadProcessedProj(
  proj = projChronic,
  adt = adtChronic,
  tsaCatalog = tsa_catalog,
  generateNewUmapSettings = list(
    "reducedDims" = "Harmony",
    "name" = "UMAP2",
    "nNeighbors" = 80,
    "force" = TRUE,
    "minDist" = 0.3),
  generateNewClusterSettings = list(
    "reducedDims" = "Harmony",
    "method" = "Seurat",
    "name" = "Clusters2",
    "resolution" = 0.5
  ),
  runPreAnnot = TRUE,
  runPreAnnotIsoComps = isoComparisonsChronic,
  prefixForGraphs = "chronic"
)

projChronic_matched <- chronicPreProcessed$archrProj
adtChronic_matched <- chronicPreProcessed$adtSeu
chronicAdtToUse <- chronicPreProcessed$adtToUse

rm(chronicPreProcessed)

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
    "force" = TRUE,
    "resolution" = 0.75
  ),
  runPreAnnot = TRUE,
  runPreAnnotIsoComps = isoComparisonsART,
  prefixForGraphs = "art"
)

projART_matched <- artPreProcessed$archrProj
adtART_matched <- artPreProcessed$adtSeu
ARTAdtToUse <- artPreProcessed$adtToUse

rm(artPreProcessed)
gc()



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
  colorScheme = scale_color_manual(values = multiHueColorPalette),
  propInLegend = TRUE)

plotDiscreteLollipop(projInvitro_matched, "inVitro_discreteAbsolute_matched",
  cluster = "manualClusterAnnot",
  graphType = "absolute")

plotDiscreteLollipop(projInvitro_matched, "inVitro_discreteHivOnly_matched",
  cluster = "manualClusterAnnot",
  graphType = "hivOnly")

###############################################################################
# chronic manual annot
###############################################################################

projChronic_matched <- assignManualAnnotation(archrProj = projChronic_matched,
  cluster = "Clusters2",
  annotFn = "manualClusterAnnotations/chronic.csv")

adtChronic_matched$manualClusterAnnot <- projChronic_matched$manualClusterAnnot
adtChronic_matched$haystackOut <- projChronic_matched$haystackOut

plotUmap(projChronic_matched,
  colorBy = "manualClusterAnnot",
  fn = "chronic_umap_labeledCluster",
  embedding = "UMAP2")
plotUmap(projChronic_matched,
  colorBy = "Clusters2",
  colorLabelBy = "Clusters2",
  fn = "chronic_umap_unannotCluster_withPlotLabels",
  embedding = "UMAP2",
  colorScheme = scale_color_manual(values = multiHueColorPalette),
  propInLegend = TRUE)
plotUmap(projChronic_matched,
  fn = "chronic_umap_labeledCluster_withPlotLabels",
  colorBy = c("Clusters2", "manualClusterAnnot"),
  colorLabelBy = "Clusters2",
  embedding = "UMAP2",
  colorScheme = scale_color_manual(values = multiHueColorPalette),
  propInLegend = TRUE)

plotDiscreteLollipop(projChronic_matched, "chronic_discreteAbsolute_matched",
  cluster = "manualClusterAnnot",
  graphType = "absolute")

plotDiscreteLollipop(projChronic_matched, "chronic_discreteHivOnly_matched",
  cluster = "manualClusterAnnot",
  graphType = "hivOnly")

###############################################################################
# art manual annot
###############################################################################
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
  colorScheme = scale_color_manual(values = multiHueColorPalette),
  propInLegend = TRUE)

plotDiscreteLollipop(projART_matched, "art_discreteAbsolute_matched",
  cluster = "manualClusterAnnot",
  gheight = 4,
  graphType = "absolute")

plotDiscreteLollipop(projART_matched, "art_discreteHivOnly_matched",
  cluster = "manualClusterAnnot",
  gheight = 4,
  donorColumn = "individual",
  graphType = "hivOnly")


plotDiscreteLollipop(projART_matched, "art_discreteAbsolute_matched_byDonor",
                     cluster = "manualClusterAnnot",
                     donorColumn = "individual",
                     gheight = 4,
                     gwidth = 14,
                     graphType = "absolute")

plotDiscreteLollipop(projART_matched, "art_discreteHivOnly_matched_byDonor",
                     cluster = "manualClusterAnnot",
                     donorColumn = "individual",
                     gheight = 4,
                     gwidth = 14,
                     graphType = "hivOnly")

# looking at the HIV+ cells in non-CD4+ t-cells
# tmp <- data.frame(
#   CD4 = adtART_matched@assays$tsa@data["A0072", ],
#   CD3 = adtART_matched@assays$tsa@data["A0034", ], 
#   CD5 = adtART_matched@assays$tsa@data["A0138", ],
#   CD14 = adtART_matched@assays$tsa@data["A0081", ],
#   CCR2 = adtART_matched@assays$tsa@data["A0242", ],
#   CD16 = adtART_matched@assays$tsa@data["A0083", ],
#   OX40 = adtART_matched@assays$tsa@data["A0158", ],
#   CD137 = adtART_matched@assays$tsa@data["A0355", ],
#   `HLA-DR` = adtART_matched@assays$tsa@data["A0159", ],
#   CD38 = adtART_matched@assays$tsa@data["A0389", ],
#   CD71 = adtART_matched@assays$tsa@data["A0394", ],
#   CD69 = adtART_matched@assays$tsa@data["A0146", ],
#   manualClusterAnnot = adtART_matched$manualClusterAnnot,
#   haystackOut = adtART_matched$haystackOut) %>%
#   pivot_longer(cols = -all_of(c("manualClusterAnnot", "haystackOut")), names_to = "marker", values_to = "expression")
# 
# ggplot(mapping = aes(x = manualClusterAnnot, y = expression)) +
#   geom_hline(yintercept = 2, color = "#cccccc", linetype = "dashed") +
#   geom_violin(data = tmp %>% filter(!haystackOut), fill = HIVNEGCOLOR, color = "#333333") +
#   geom_point(data = tmp %>% filter(haystackOut), shape = 21, fill = HIVPOSCOLOR, color = "#000000") +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   facet_wrap(~ marker, ncol = 1, scales = "free_y", strip.position = "right") +
#   labs(x = "Cluster", y = "Expression")


###############################################################################
# iv differential analysis
###############################################################################
invitro_tcells_check <- grepl("^CD4", adtInvitro_matched$manualClusterAnnot)
invitro_tcells <- adtInvitro_matched$haystackOut[invitro_tcells_check]
invitro_markers <- findHIVDifferentialMarkers(seu = adtInvitro_matched,
  cellsPos = names(invitro_tcells[invitro_tcells]),
  cellsNeg = names(invitro_tcells[!invitro_tcells]),
  featuresToUse = invitroAdtToUse,
  tsa_catalog = tsa_catalog)

makeDifferentialLollipop(invitro_markers, "invitro_lollipop")
exportTsv(invitro_markers)

invitro_activeLaterCells <- adtInvitro_matched$haystackOut[grepl("CD4 Activated", adtInvitro_matched$manualClusterAnnot) & invitro_tcells_check]
invitro_activeLate_markers <- findHIVDifferentialMarkers(
  seu = adtInvitro_matched,
  tsa_catalog = tsa_catalog,
  cellsPos = names(invitro_activeLaterCells[invitro_activeLaterCells]),
  cellsNeg = names(invitro_activeLaterCells[!invitro_activeLaterCells]),
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
  cellsPos = names(invitro_earlyCells[invitro_earlyCells]),
  cellsNeg = names(invitro_earlyCells[!invitro_earlyCells]),
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
    minReplicates = 3,
    minCells = 333,
    force = TRUE)
  
  projInvitro_matched <- addReproduciblePeakSet(
    ArchRProj = projInvitro_matched,
    groupBy = "condensedHivCluster", 
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

invitro_markerTest_early <- getMarkerFeatures(
  ArchRProj = projInvitro_matched_gcPeak, 
  useMatrix = "PeakMatrix",
  groupBy = "condensedHivCluster",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = c("Naive_earlyMemory_TRUE", "Tcm_TRUE"),
  bgdGroups = c("Naive_earlyMemory_FALSE", "Tcm_FALSE"),
  maxCells = 1000,
  verbose = FALSE
)

plotVolcanoFromGetMarkerFeatures(invitro_markerTest_early,
  fn = "invitro_peaks_volcano_earlyHIV")

# chromVAR analysis
invitro_chromvar_markerTest <- getMarkerFeatures(
  ArchRProj = projInvitro_matched_gcPeak,
  useMatrix = "MotifMatrix",
  groupBy = "condensedHivCluster",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Activated_TRUE",
  bgdGroups = "Activated_FALSE",
  maxCells = 1000,
  verbose = FALSE,
  useSeqnames = "z"
)

plotVolcanoFromGetMarkerFeatures(invitro_chromvar_markerTest, chromVARmode = TRUE,
  fn = "invitro_volcano_motifs_chromVAR")

plotMotifDot(invitro_chromvar_markerTest, "invitro_chromVAR_motifsUp_activatedHIVneg", direction = "negative")
plotMotifDot(invitro_chromvar_markerTest, "invitro_chromVAR_motifsUp_activatedHIVpos", direction = "positive")

invitro_chromvar_markerTest_early <- getMarkerFeatures(
  ArchRProj = projInvitro_matched_gcPeak,
  useMatrix = "MotifMatrix",
  groupBy = "condensedHivCluster",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Tcm_TRUE",
  bgdGroups = "Tcm_FALSE",
  verbose = FALSE,
  useSeqnames = "z"
)

plotVolcanoFromGetMarkerFeatures(invitro_chromvar_markerTest_early, chromVARmode = TRUE,
  fn = "invitro_volcano_motifs_early_chromVAR")

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
# chronic differential analysis
###############################################################################
chronic_tcells_check <- !(grepl("(B|APC|CD8)", adtChronic_matched$manualClusterAnnot))
chronic_tcells <- adtChronic_matched$haystackOut[chronic_tcells_check]
chronic_markers <- findHIVDifferentialMarkers(seu = adtChronic_matched,
  cellsPos = names(chronic_tcells[chronic_tcells]),
  cellsNeg = names(chronic_tcells[!chronic_tcells]),
  featuresToUse = chronicAdtToUse,
  findMarkerMethod = "DESeq2",
  tsa_catalog = tsa_catalog)

plotVlnEnhanced(adtChronic_matched, cells = names(chronic_tcells),
  feats = chronic_markers$gene,
  separator = chronic_markers$Status,
  titles = chronic_markers$cleanName,
  fn = "chronic_tcells_vln")

makeDifferentialLollipop(chronic_markers, "chronic_lollipop")
exportTsv(chronic_markers)

chronic_tfh_cells <- adtChronic_matched$haystackOut[grepl("(CD4 Tfh)", adtChronic_matched$manualClusterAnnot)]
chronic_tfh_markers <- findHIVDifferentialMarkers(
  seu = adtChronic_matched,
  tsa_catalog = tsa_catalog,
  cellsPos = names(chronic_tfh_cells[chronic_tfh_cells]),
  cellsNeg = names(chronic_tfh_cells[!chronic_tfh_cells]),
  findMarkerMethod = "DESeq2",
  featuresToUse = chronicAdtToUse)

makeDifferentialLollipop(chronic_tfh_markers, "chronic_lollipop_tfh")
exportTsv(chronic_tfh_markers)

plotVlnEnhanced(adtChronic_matched, cells = names(chronic_tfh_cells),
  feats = chronic_tfh_markers$gene,
  separator = chronic_tfh_markers$Status,
  titles = chronic_tfh_markers$cleanName,
  fn = "chronic_tfh_vln")

chronicGcPeakFn <- paste0(getOutputDirectory(projChronic_matched), "_gcPeakCalled")
if (!dir.exists(chronicGcPeakFn)) {
  chronic_renamingDF <- getCellColData(projChronic_matched, select = c("haystackOut", "manualClusterAnnot")) %>%
    as.data.frame(.) %>%
    mutate(cellNames = projChronic_matched$cellNames) %>%
    mutate(condensedCluster = case_when(
      grepl("^CD4", manualClusterAnnot) ~ "CD4",
      TRUE ~ "Other")) %>%
    mutate(condensedHIVCluster = paste0(condensedCluster, "_", haystackOut))
  
  projChronic_matched <- addCellColData(ArchRProj = projChronic_matched,
    data = chronic_renamingDF$condensedHIVCluster,
    cells = chronic_renamingDF$cellNames,
    name = "condensedHivCluster",
    force = TRUE)
  
  projChronic_matched <- addGroupCoverages(ArchRProj = projChronic_matched,
    groupBy = "condensedHivCluster",
    minReplicates = 3,
    minCells = 15,
    sampleRatio = 0.95,
    force = TRUE)
  
  projChronic_matched <- addReproduciblePeakSet(
    ArchRProj = projChronic_matched,
    groupBy = "condensedHivCluster", 
    pathToMacs2 = "/Users/vince/anaconda3/envs/asapseq/bin/macs2",
    force = TRUE)
  
  projChronic_matched <- addPeakMatrix(projChronic_matched, force = TRUE)
  
  if("Motif" %ni% names(projChronic_matched@peakAnnotation)) {
    projChronic_matched <- addMotifAnnotations(ArchRProj = projChronic_matched,
      motifSet = "cisbp", name = "Motif", force = TRUE)
  }
  
  projChronic_matched <- addBgdPeaks(projChronic_matched, force = TRUE)
  
  projChronic_matched <- addDeviationsMatrix(
    ArchRProj = projChronic_matched, 
    peakAnnotation = "Motif",
    force = TRUE
  )
  
  projChronic_matched_gcPeak <- saveArchRProject(projChronic_matched,
    outputDirectory = chronicGcPeakFn, load = TRUE)
} else {
  projChronic_matched_gcPeak <- loadArchRProject(chronicGcPeakFn, force = TRUE, showLogo = FALSE)
}


chronic_markerTest <- getMarkerFeatures(
  ArchRProj = projChronic_matched_gcPeak, 
  useMatrix = "PeakMatrix",
  groupBy = "condensedHivCluster",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CD4_TRUE",
  bgdGroups = "CD4_FALSE",
  maxCells = 1000,
  verbose = FALSE
)
plotVolcanoFromGetMarkerFeatures(chronic_markerTest,
  fn = "chronic_peaks_volcano")

chronic_markerTest_chromvar <- getMarkerFeatures(
  ArchRProj = projChronic_matched_gcPeak,
  useMatrix = "MotifMatrix",
  groupBy = "condensedHivCluster",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CD4_TRUE",
  bgdGroups = "CD4_FALSE",
  maxCells = 1000,
  verbose = TRUE,
  useSeqnames = "z"
)

plotVolcanoFromGetMarkerFeatures(chronic_markerTest_chromvar, chromVARmode = TRUE,
  fn = "chronic_volcano_motifs_chromVAR")

plotMotifDot(chronic_markerTest_chromvar, "chronic_chromVAR_motifsUp_HIVneg", direction = "negative")
plotMotifDot(chronic_markerTest_chromvar, "chronic_chromVAR_motifsUp_HIVpos", direction = "positive")

###############################################################################
# chronic model analysis
###############################################################################
chronicModelData <- getAssayData(adtChronic_matched, names(chronic_tcells), chronicAdtToUse, splitRatio = 0.7)
chronicLR <- calcLogisticalRegression(chronicModelData[["train"]])
chronicNB <- calcNaiveBayes(chronicModelData[["train"]])
chronicRF <- calcRandomForest(chronicModelData[["train"]],
  sampleProps = list(c(NULL, NULL), c(5, 1), c(10, 1)))

chronicRFPred <- lapply(chronicRF, function(x) {
  pred <- getClassifierPrediction(predict(x, chronicModelData[["test"]], type = "prob")[, 2],
    chronicModelData[["test"]]$haystackOut)
  
  return(pred)
})
names(chronicRFPred) <- names(chronicRF)

chronicModelsPred <- list(
  "logistic" = getClassifierPrediction(predict(chronicLR, chronicModelData[["test"]], type = 'response'),
    chronicModelData[["test"]]$haystackOut),
  "naive bayes" = getClassifierPrediction(predict(chronicNB, chronicModelData[["test"]], type = "prob")[, 2],
    chronicModelData[["test"]]$haystackOut)
)

chronicModelsPred <- base::append(chronicModelsPred, chronicRFPred)
chronicModelsROC <- lapply(chronicModelsPred, getClassifierPerformance)
chronicModelsAUC <- lapply(chronicModelsPred, getClassifierPerformance, perfMode = "auc")

plotGgRoc(chronicModelsROC, chronicModelsAUC, "chronic_models")
plotLogitRegressionVolcano(chronicLR, tsa_catalog, chronicAdtToUse, "chronic_logitRegression_volcano")

randomForest::importance(chronicRF$`RF (all)`, type = "1") %>% 
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
art_markers <- findHIVDifferentialMarkers(seu = adtART_matched,
  cellsPos = names(art_tcells[art_tcells]),
  cellsNeg = names(art_tcells[!art_tcells]),
  featuresToUse = ARTAdtToUse,
  findMarkerMethod = "DESeq2",
  tsa_catalog = tsa_catalog)

makeDifferentialLollipop(art_markers, "art_lollipop")
exportTsv(art_markers)

plotVlnEnhanced(adtART_matched, cells = names(art_tcells),
  feats = art_markers$gene,
  separator = art_markers$Status,
  titles = art_markers$cleanName,
  fn = "art_tcells_vln")

tmp <- data.frame(CD4 = adtART_matched@assays$tsa@data["A0072", ],
  CD3 = adtART_matched@assays$tsa@data["A0034", ], 
  CD5 = adtART_matched@assays$tsa@data["A0138", ],
  CD14 = adtART_matched@assays$tsa@data["A0081", ],
  manualClusterAnnot = adtART_matched$manualClusterAnnot,
  haystackOut = adtART_matched$haystackOut) %>%
  pivot_longer(cols = c("CD4", "CD3", "CD5", "CD14"), names_to = "marker", values_to = "expression")

(ggplot(mapping = aes(x = manualClusterAnnot, y = expression)) +
    geom_hline(yintercept = 2, color = "#cccccc") +
    geom_violin(data = tmp %>% filter(!haystackOut), fill = HIVNEGCOLOR, color = "#333333") +
    geom_point(data = tmp %>% filter(haystackOut), shape = 21, fill = HIVPOSCOLOR, color = "#000000") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~ marker, ncol = 1, scales = "free_y", strip.position = "right") +
    labs(x = "Cluster", y = "Expression")) %>%
  ggsave(filename = "outs/png/art_hivVln.png", plot = ., dpi = "retina", width = 7, height = 5)


artGcPeakFn <- paste0(getOutputDirectory(projART_matched), "_gcPeakCalled")
if (!dir.exists(artGcPeakFn)) {
  art_renamingDF <- getCellColData(projART_matched, select = c("haystackOut", "manualClusterAnnot")) %>%
    as.data.frame(.) %>%
    mutate(cellNames = projART_matched$cellNames) %>%
    mutate(condensedCluster = case_when(
      grepl("CD4", manualClusterAnnot) ~ "T cell",
      TRUE ~ "Other")) %>%
    mutate(condensedHIVCluster = paste0(condensedCluster, "_", haystackOut))
  
  projART_matched <- addCellColData(ArchRProj = projART_matched,
    data = art_renamingDF$condensedHIVCluster,
    cells = art_renamingDF$cellNames,
    name = "condensedHivCluster",
    force = TRUE)
  
  projART_matched <- addGroupCoverages(ArchRProj = projART_matched,
    groupBy = "condensedHivCluster",
    minReplicates = 2,
    minCells = 18,
    sampleRatio = 0.95,
    force = TRUE)
  
  projART_matched <- addReproduciblePeakSet(
    ArchRProj = projART_matched,
    groupBy = "condensedHivCluster", 
    pathToMacs2 = "/Users/vince/anaconda3/envs/asapseq/bin/macs2",
    force = TRUE)
  
  projART_matched <- addPeakMatrix(projART_matched, force = TRUE)
  
  if ("Motif" %ni% names(projART_matched@peakAnnotation)) {
    projART_matched <- addMotifAnnotations(ArchRProj = projART_matched, motifSet = "cisbp", name = "Motif")
  }
  
  projART_matched <- addBgdPeaks(projART_matched, force = TRUE)
  
  projART_matched <- addDeviationsMatrix(
    ArchRProj = projART_matched, 
    peakAnnotation = "Motif",
    force = TRUE
  )
  
  if ("EncodeTFBS" %ni% names(projART_matched@peakAnnotation)) {
    projART_matched <- addArchRAnnotations(ArchRProj = projART_matched, collection = "EncodeTFBS")
  }
  
  projART_matched <- addDeviationsMatrix(
    ArchRProj = projART_matched, 
    peakAnnotation = "EncodeTFBS",
    force = TRUE
  )
  
  projART_matched_gcPeak <- saveArchRProject(projART_matched, outputDirectory = artGcPeakFn, load = TRUE)
} else {
  projART_matched_gcPeak <- loadArchRProject(artGcPeakFn)
}

art_markerTest <- getMarkerFeatures(
  ArchRProj = projART_matched_gcPeak, 
  useMatrix = "PeakMatrix",
  groupBy = "condensedHivCluster",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "T cell_TRUE",
  bgdGroups = "T cell_FALSE",
  maxCells = 1000,
  verbose = FALSE
)

plotVolcanoFromGetMarkerFeatures(art_markerTest, chromVARmode = FALSE, fn = "art_volcano")

art_markerTest_deviations <- getMarkerFeatures(
  ArchRProj = projART_matched_gcPeak, 
  useMatrix = "MotifMatrix",
  groupBy = "condensedHivCluster",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "T cell_TRUE",
  bgdGroups = "T cell_FALSE",
  maxCells = 500,
  verbose = FALSE,
  useSeqnames = "z",
)

plotVolcanoFromGetMarkerFeatures(art_markerTest_deviations, chromVARmode = TRUE,
  fn = "art_volcano_motifs_chromVAR")
plotMotifDot(art_markerTest_deviations, "art_chromVAR_motifsUp_HIVneg", direction = "negative")
plotMotifDot(art_markerTest_deviations, "art_chromVAR_motifsUp_HIVpos", direction = "positive")

art_markerTest_encode_deviations <- getMarkerFeatures(
  ArchRProj = projART_matched_gcPeak, 
  useMatrix = "EncodeTFBSMatrix",
  groupBy = "condensedHivCluster",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "T cell_TRUE",
  bgdGroups = "T cell_FALSE",
  maxCells = 500,
  verbose = TRUE,
  useSeqnames = "z"
)

plotVolcanoFromGetMarkerFeatures(art_markerTest_encode_deviations, chromVARmode = TRUE,
  fn = "art_volcano_motifs_chromVAR_encode")
plotMotifDot(art_markerTest_encode_deviations, "art_chromVAR_encodeUp_HIVneg", direction = "negative")
plotMotifDot(art_markerTest_encode_deviations, "art_chromVAR_encodeUp_HIVpos", direction = "positive")




###############################################################################
# qc and frags
###############################################################################
plotUpsetCBC(atacCBC = projInVitro$cellNames,
  hivCBC = unique(haystackInVitro$allViralFrags$newCbc),
  adtCBC = Cells(adtInVitro),
  fn = "invitro_upset_cbc")

plotUpsetCBC(atacCBC = projChronic$cellNames,
  hivCBC = unique(haystackChronic$allViralFrags$newCbc),
  adtCBC = Cells(adtChronic),
  separateByIndividual = TRUE,
  fn = "chronic_upset_cbc")

plotUpsetCBC(atacCBC = projART$cellNames,
  hivCBC = unique(haystackART$allViralFrags$newCbc),
  adtCBC = Cells(adtART),
  separateByIndividual = TRUE,
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

chronicFrags_matched <- haystackChronic$filteredviralFrags[haystackChronic$filteredviralFrags$newCbc %in% projChronic_matched$cellNames, ] %>%
  mutate(individual = ifelse(grepl("C01", sample), "C01", "C02"))

plotFragMultiGraph(frags = chronicFrags_matched, coverage = rep(0, 9719), geneAnnots = hxb2GeneAnnots, separateByIndividual = TRUE, fn = "chronic_frags", gheight = 6)
calculateLTRVersusInternal(chronicFrags_matched, 633, 9085, 9718)

artFrags_matched <- haystackART$filteredviralFrags[haystackART$filteredviralFrags$newCbc %in% projART_matched$cellNames, ]

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

artConsensusAnnots <- read.csv("viralGenomeAnnotations/GeneCutterParser_A08_A01_B45.tsv", sep = "\t") %>%
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

makeFancyUpsetPlotHIV(
  seu = adtInvitro_matched,
  cellsOfInterest = names(invitro_activeLaterCells),
  featuresOfInterest = c("A0141", "A0870", "A0367", "A0147", "A0390"),
  featuresOfInterestNames = c("CCR5", "SLAM", "CD2", "CD62L", "CD127"),
  metadata = ifelse(invitro_activeLaterCells, "HIV+", "HIV-"),
  thresholds = c(1.25, 1.25, 1.25, 2.75, 1.5),
  fn = "invitro_fancyUpset_activatedLaterCells"
)

tmp <- calculateCoverage(inVitroFrags_matched, rep(0, 9726))

tmp2 <- calculateCoverage(inVitroFrags_matched %>%
    group_by(readname, newCbc) %>%
    filter(n() == 2) %>%
    arrange(min(startBp, endBp), .by_group = TRUE) %>%
    summarise(inferredStart = first(endBp) + 1, inferredEnd = last(startBp) -1) %>%
    as.data.frame(.),
  coverage = rep(0, 9726), startCol = "inferredStart", endCol = "inferredEnd")

tmp %>%
  mutate(y = y + tmp2$y) %>%
  filter(x < 800) %>%
  ggplot(aes(x = x, y=y)) +
  annotate("rect",
    xmin = c(40, 456-2), xmax = c(200, 456+140), ymin = 0, ymax = 0.3,
    alpha = .1, fill = "blue") +
  annotate("rect",
    xmin = c(315), xmax = c(408), ymin = 0, ymax = 0.3,
    alpha = .1, fill = "red") +
  geom_vline(xintercept = 456, color = "red") +
  geom_line() +
  theme_classic()


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
  

