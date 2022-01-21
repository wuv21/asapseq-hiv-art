library(glue)

exportTsv <- function(
  df,
  fn = deparse(substitute(df)),
  motifMode = FALSE,
  chromVarMode = FALSE) {
  
  if (motifMode & chromVarMode) {
    df <- data.frame(
      y = log10(assays(df)$FDR[, 1]) * -1,
      meandiff = assays(df)$MeanDiff[, 1],
      motif = elementMetadata(df)$name)

  } else if (motifMode) {
    df <- data.frame(
      motif = rownames(df),
      y = assay(df)[, 1])
  }
  
  write.table(df, file = glue("outs/tsv/{fn}.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE)
}