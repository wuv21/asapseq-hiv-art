library(rentrez)

getFeatureTable <- function(id) {
  res <- entrez_fetch(db = "nuccore", id = id, rettype = "ft")
  
  features <- str_split(res, "\n\t\t\tgene\t")[[1]]
  features2 <- str_match(features, "(\\w+)\n<?(\\d+)\t(\\d+)\t(CDS|misc_feature)\n(\\d*)\t(\\d*)")
  
  df <- data.frame(annotation = features2[-1, 2],
    startPos = as.numeric(features2[-1, 3]) - 1,
    endPos = as.numeric(features2[-1, 4]) - 1,
    startPos2 = as.numeric(features2[-1, 6]),
    endPos2 = as.numeric(features2[-1, 7])) %>%
    mutate(startPos2 = ifelse(is.na(startPos2), NA, startPos2 - 1),
      endPos2 = ifelse(is.na(endPos2), NA, endPos2 - 1))
  
  dfExtra <- df[, c("annotation", "startPos2", "endPos2")] %>%
    dplyr::filter(!is.na(startPos2)) %>%
    dplyr::rename(startPos = startPos2, endPos = endPos2)
  
  df <- df %>%
    dplyr::select(-startPos2, -endPos2) %>%
    dplyr::bind_rows(dfExtra)
  
  return(df)
}

downloadTables <- function(ids) {
  allTables <- list()
  for (id in ids) {
    if (grepl("A08", id)){
      fn <- glue("viralGenomeAnnotations/{id}.csv")
      
      if (file.exists(fn)) {
        allTables[[id]] <- read.csv(fn)
      } else {
        idClean <- gsub("chrA08", "", id)
        featureDf <- getFeatureTable(idClean)
        featureDf$genome <- rep(id, nrow(featureDf))
        
        write.table(x = featureDf, file = fn, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")
        allTables[[id]] <- featureDf
      }
    }
  }
  
  return(allTables)
}