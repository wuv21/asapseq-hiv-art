library(tidyverse)
library(glue)

tsa_catalog <- read.csv("tsa_catalog/tsa_catalog_raw.csv", header = TRUE)

tsa_catalog <- tsa_catalog %>%
  mutate(isCtrl = grepl("Isotype Ctrl", Description) | grepl("isotype Ctrl", Description),
    cleanDescriptionFull = gsub("anti-(\\w|\\/)+ ", "", Description),
    nonCDName = str_match(Description, "\\((.*)\\)")[, 2],
    cdValue = str_match(Description, "(CD\\d+\\w*)")[, 2],
    cleanName = ifelse(is.na(cdValue),
      cleanDescriptionFull, ifelse(is.na(nonCDName), cdValue, glue("{cdValue} ({nonCDName})")))) 

saveRDS(tsa_catalog, file = "rds/tsa_catalog.rds")
