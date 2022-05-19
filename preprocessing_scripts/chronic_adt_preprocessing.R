suppressMessages(library(Seurat))
suppressMessages(library(DropletUtils))
suppressMessages(library(tidyverse))
suppressMessages(library(glue))
source("adt_preprocessing_general.R")

set.seed(21)

samples <- c("C01_1", "C01_2", "C01_3", "C02_1", "C02_2", "C02_3")

tsa_catalog <- readRDS("../rds/tsa_catalog.rds")

seuListFn <-"../rds/chronic_seuratList.rds"
if (!file.exists(seuListFn)) {
  seuList <- createSeuList(samples = samples, adtDir = "../data/adt", filename = seuListFn)
} else {
  seuList <- readRDS(file = seuListFn)
}

adtChronic <- merge(x = seuList[[1]],
  y = seuList[2:length(seuList)],
  add.cell.ids = paste0(samples, "#"),
  merge.data = TRUE)

adtChronic <- NormalizeData(adtChronic, normalization.method = "CLR", scale.factor = 10000, margin = 2)
adtChronic <- ScaleData(adtChronic) %>%
  RunPCA(features = rownames(adtChronic@assays$tsa))

adtChronic$individual <- str_match(adtChronic$orig.ident, "C\\d+")[,1]

adtChronic <- adtChronic %>%
  harmony::RunHarmony(group.by.vars = "individual", plot_convergence = TRUE, assay.use = "tsa")

adtChronic <- RenameCells(adtChronic, new.names = gsub("#_", "#", paste0(Cells(adtChronic), "-1")))

isoControls <- tsa_catalog[tsa_catalog$isCtrl, ]
nonIsoControls <- tsa_catalog[!tsa_catalog$isCtrl, ]

isoComparisonsChronic <- lapply(nonIsoControls$DNA_ID, function(nonIsoId) {
  comps <- sapply(isoControls$DNA_ID, function(isoID) {
    nonIsoCounts <- adtChronic@assays$tsa@data[nonIsoId, ]
    isoCounts <- adtChronic@assays$tsa@data[isoID, ]
    
    return(wilcox.test(nonIsoCounts, isoCounts, alternative = "greater")$p.value)
  })
  
  names(comps) <- isoControls$DNA_ID
  return(comps)
})

isoComparisonsChronic <- data.frame(bind_rows(isoComparisonsChronic))
rownames(isoComparisonsChronic) <- nonIsoControls$DNA_ID

save(isoComparisonsChronic, adtChronic, file = "../rds/chronic_seuratMerged.RData")
