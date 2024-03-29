suppressMessages(library(Seurat))
suppressMessages(library(DropletUtils))
suppressMessages(library(tidyverse))
suppressMessages(library(glue))
source("adt_preprocessing_general.R")

set.seed(21)

samples <- c(
  "A01_1", "A01_2", "A01_3", "A01_4",
  "A01_pre_1", "A01_pre_2", "A01_pre_3",
  "A08_1", "A08_2", "A08_4",
  "A08_pre_1", "A08_pre_2", "A08_pre_3", "A08_pre_4",
  "BEAT045_A2", "BEAT045_A3", "BEAT045_B1", "BEAT045_B2",
  "A09_pre_1", "A09_pre_2", "A09_pre_3", "A09_pre_4", "A09_pre_5", "A09_pre_6",
  "A09_post_2", "A09_post_3", "A09_post_4", "A09_post_5", "A09_post_6")

tsa_catalog <- readRDS("../rds/tsa_catalog.rds")

seuListFn <-"../rds/ART_seuratList.rds"
if (!file.exists(seuListFn)) {
  seuList <- createSeuList(samples = samples, adtDir = "../data/adt", filename = seuListFn, emptyDropsLower = 250)
} else {
  seuList <- readRDS(file = seuListFn)
}

adtART <- merge(x = seuList[[1]],
  y = seuList[2:length(seuList)],
  add.cell.ids = paste0(samples, "#"),
  merge.data = TRUE)

adtART <- RenameCells(adtART, new.names = gsub("BEAT045", "B45", Cells(adtART)))
adtART$orig.ident <- gsub("BEAT045", "B45", adtART$orig.ident)

# change a01 and a08 original samples to a01_post and a08_post respectively
adtART <- RenameCells(adtART, new.names = stringr::str_replace(Cells(adtART), "(A0\\d)_(\\d)", "\\1_post_\\2"))
adtART$orig.ident <- stringr::str_replace(adtART$orig.ident, "(A0\\d)_(\\d)", "\\1_post_\\2")

adtART <- NormalizeData(adtART, normalization.method = "CLR", scale.factor = 10000, margin = 2)
adtART <- ScaleData(adtART) %>%
  RunPCA(features = rownames(adtART@assays$tsa))

adtART$individual <- stringr::str_match(adtART$orig.ident, "\\w\\d{2}(_pre|_post)*")[, 1]

adtART <- adtART %>%
  harmony::RunHarmony(group.by.vars = "orig.ident", plot_convergence = TRUE, assay.use = "tsa")

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

isoComparisonsART <- data.frame(bind_rows(isoComparisonsART))
rownames(isoComparisonsART) <- nonIsoControls$DNA_ID

save(isoComparisonsART, adtART, file = "../rds/ART_seuratMerged.RData")
