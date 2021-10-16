suppressMessages(library(Seurat))
suppressMessages(library(DropletUtils))
suppressMessages(library(tidyverse))
suppressMessages(library(glue))
source("adt_preprocessing_general.R")

samples <- c("A01_1", "A01_2", "A01_3", "A01_4", "A08_1", "A08_2", "A08_4")

tsa_catalog <- readRDS("../rds/tsa_catalog.rds")

seuListFn <-"../rds/ART_seuratList.rds"
if (!file.exists(seuListFn)) {
  seuList <- createSeuList(samples = samples, adtDir = "../data/adt", filename = seuListFn)
} else {
  seuList <- readRDS(file = seuListFn)
}

adtART <- merge(x = seuList[[1]],
  y = seuList[2:length(seuList)],
  add.cell.ids = paste0(samples, "#"),
  merge.data = TRUE)

adtART <- NormalizeData(adtART, normalization.method = "CLR", scale.factor = 10000, margin = 2)
adtART <- ScaleData(adtART) %>%
  RunPCA(features = rownames(adtART@assays$tsa))

adtART$individual <- str_match(adtART$orig.ident, "A\\d+")[,1]

adtART <- adtART %>%
  harmony::RunHarmony(group.by.vars = "individual", plot_convergence = TRUE)

adtART <- RenameCells(adtART, new.names = gsub("#_", "#", paste0(Cells(adtART), "-1")))

isoControls <- tsa_catalog[tsa_catalog$isCtrl, ]
nonIsoControls <- tsa_catalog[!tsa_catalog$isCtrl, ]

isoComparisonsART <- lapply(nonIsoControls$DNA_ID, function(nonIsoId) {
  comps <- sapply(isoControls$DNA_ID, function(isoID) {
    nonIsoCounts <- adtART@assays$tsa@data[nonIsoId, ]
    isoCounts <- adtART@assays$tsa@data[isoID, ]
    
    return(wilcox.test(nonIsoCounts, isoCounts, alternative = "greater")$p.value)
  })
  
  names(comps) <- isoControls$DNA_ID
  return(comps)
})

save(isoComparisonsART, adtART, file = "../rds/ART_seuratMerged.RData")
