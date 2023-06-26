
#' Return function to interpolate a continuous POMA color palette
#'
#' @param palette Character name of palette in poma_palettes
#'
#' @export
poma_pal_c <- function(palette = "nature") {
  poma_palettes <- list(
    `nature`      = c("#E64B35FF", "#4DBBD5FF", "#00A087FF"),
    `simpsons`    = c("#FED439FF", "#709AE1FF", "#8A9197FF"),
    `futurama`    = c("#FF6F00FF", "#C71000FF", "#008EA0FF")
  )
  
  pal <- poma_palettes[[palette]]
}

#' Return function to interpolate a discrete POMA color palette
#'
#' @param palette Character name of palette in poma_palettes
#'
#' @export
poma_pal_d <- function(palette = "nature") {
  poma_palettes <- list(
    `nature`      = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF"),
    `simpsons`    = c("#FED439FF", "#709AE1FF", "#8A9197FF", "#D2AF81FF", "#FD7446FF", "#D5E4A2FF", "#197EC0FF", "#F05C3BFF", "#46732EFF", "#71D0F5FF"),
    `futurama`    = c("#FF6F00FF", "#C71000FF", "#008EA0FF", "#8A4198FF", "#5A9599FF", "#FF6348FF", "#84D7E1FF", "#FF95A8FF", "#3D3B25FF", "#ADE2D0FF")
  )

  pal <- poma_palettes[[palette]]
}

#' Color scale constructor for continuous POMA color palettes
#'
#' @param palette Character name of palette in poma_palettes
#' 
# scale_color_poma_c <- function(palette = "nature") {
#   ggplot2::scale_colour_gradientn(colours = poma_pal_c(palette = palette))
# }

#' Color scale constructor for discrete POMA colors
#'
#' @param palette Character name of palette in poma_palettes
#'
# scale_color_poma_d <- function(palette = "nature") {
#   ggplot2::scale_color_manual(values = poma_pal_d(palette = palette))
# }

#' Fill scale constructor for continuous POMA color palettes
#'
#' @param palette Character name of palette in poma_palettes
#'
# scale_fill_poma_c <- function(palette = "nature") {
#   ggplot2::scale_fill_gradientn(colours = poma_pal_c(palette = palette))
# }

#' Fill scale constructor for discrete POMA color palettes
#'
#' @param palette Character name of palette in poma_palettes
#'
# scale_fill_poma_d <- function(palette = "nature") {
#   ggplot2::scale_fill_manual(values = poma_pal_d(palette = palette))
# }

#' Color scale constructor for discrete `viridis` "plasma" palette
#'
#' @export
scale_color_poma_d <- function() {
  ggplot2::scale_color_viridis_d(option = "plasma", end = 0.8)
}

#' Fill scale constructor for discrete `viridis` "plasma" palette
#'
#' @export
scale_fill_poma_d <- function() {
  ggplot2::scale_fill_viridis_d(option = "plasma", end = 0.8)
}

#' Color scale constructor for continuous `viridis` "plasma" palette
#'
#' @export
scale_color_poma_c <- function() {
  ggplot2::scale_color_viridis_c(option = "plasma", end = 0.8)
}

#' Fill scale constructor for continuous `viridis` "plasma" palette
#'
#' @export
scale_fill_poma_c <- function() {
  ggplot2::scale_fill_viridis_c(option = "plasma", end = 0.8)
}

