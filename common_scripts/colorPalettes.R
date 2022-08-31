library(Polychrome)

multiHueColorPalette <- c(
  "#016eb7",
  "#48b200",
  "#cc04d3",
  "#91b846",
  "#875bff",
  "#ee9300",
  "#5600a2",
  "#e69b44",
  "#027ef6",
  "#ff8d3d",
  "#322c8e",
  "#ff4549",
  "#4fb5e8",
  "#ff326b",
  "#004821",
  "#ff178c",
  "#acaf73",
  "#c895fa",
  "#4d3d00",
  "#004b86",
  "#ff7d6c",
  "#363472",
  "#751d00")

nColorPalette <- function(n, seed = 9, initcolors = c("#ff0000", "#00ff00", "#0000ff")) {
  set.seed(seed)
  p <- Polychrome::createPalette(n, seedcolors = initcolors)
  
  return(unname(p))
}