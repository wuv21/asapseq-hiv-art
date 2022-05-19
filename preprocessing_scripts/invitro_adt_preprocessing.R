suppressMessages(library(Seurat))
suppressMessages(library(DropletUtils))
suppressMessages(library(tidyverse))
suppressMessages(library(glue))
source("adt_preprocessing_general.R")

set.seed(21)

samples <- c("ND497_B")

tsa_catalog <- readRDS("../rds/tsa_catalog.rds")

isoControls <- tsa_catalog[tsa_catalog$isCtrl, ]
nonIsoControls <- tsa_catalog[!tsa_catalog$isCtrl, ]

seuListFn <-"../rds/invitro_seuratList.rds"
if (!file.exists(seuListFn)) {
  seuList <- createSeuList(samples = samples, adtDir = "../data/adt", filename = seuListFn)
  # seuList <- createSeuListWithDsb(samples = samples, adtDir = "../data/adt",
  #                                 filename = seuListFn,
  #                                 isotypeControlNames = tsa_catalog$DNA_ID)
} else {
  seuList <- readRDS(file = seuListFn)
}

adtInVitro <- seuList[[1]]
newCellNames <- paste0("ND497_B#", gsub("#_", "#", paste0(Cells(adtInVitro), "-1")))
adtInVitro <- RenameCells(adtInVitro,
  new.names = newCellNames)

isoComparisonsInvitro <- lapply(nonIsoControls$DNA_ID, function(nonIsoId) {
  comps <- sapply(isoControls$DNA_ID, function(isoID) {
    nonIsoCounts <- adtInVitro@assays$tsa@data[nonIsoId, ]
    isoCounts <- adtInVitro@assays$tsa@data[isoID, ]
    
    return(wilcox.test(nonIsoCounts, isoCounts, alternative = "greater")$p.value)
  })
  
  names(comps) <- isoControls$DNA_ID
  return(comps)
})

isoComparisonsInvitro <- data.frame(bind_rows(isoComparisonsInvitro))
rownames(isoComparisonsInvitro) <- nonIsoControls$DNA_ID

save(isoComparisonsInvitro, adtInVitro, file = "../rds/invitro_seuratMerged.RData")
