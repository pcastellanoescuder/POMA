
#' A ggplot theme which allow custom yet consistent styling of plots in the
#' POMA package and web app.
#'
#' @param base_size (integer) Base point size
#' @param axistitle (string) Axis titles. Options include "none" or
#' any combination of "X", "Y", "x" and "y".
#' @param axistext (string) Axis text labels for values or groups.
#' Options include "none" or any combination of "X", "Y", "x" and "y".
#' @param legend_position Character. Legend position. See `ggplot2` documentation.
#' @param legend_title Logical. Include legend title.
#' @param axis_x_rotate Logical. Rotate x-axis 45 degrees.
#' @param margin (numeric) Should a margin of x be added to the plot?
#' Defaults to 0 (no margin by default).
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(diamonds, aes(cut)) + geom_bar() + theme_poma()
#' }
#' 
#' @export
theme_poma <- function(base_size = 15,
                       axistitle = "xy",
                       axistext = "xy",
                       legend_position = "bottom",
                       legend_title = TRUE,
                       axis_x_rotate = FALSE,
                       margin = 2) {
  
  if(!is.character(axistitle)) stop('axistitle must be a character: "none" or any combination of "X", "Y", "x" and "y"')
  if(!is.character(axistext)) stop('axistext must be a character: "none" or any combination of "X", "Y", "x" and "y"')
  if(!is.numeric(margin)) stop('margin must be a numeric value')
  
  fontfamily_slab <- ""
  fontfamily_mono <- ""
  
  base_col <- "black"
  light_col <- "grey20"
  
  out <-
    ggplot2::theme_bw(base_size = base_size, base_family = fontfamily_slab) +
    ggplot2::theme(
      text = ggplot2::element_text(
        color = base_col
      ),
      line = ggplot2::element_line(
        color = light_col
      ),
      rect = ggplot2::element_rect(
        color = light_col,
        fill = "transparent"
      ),
      plot.title = ggtext::element_textbox_simple(
        family = fontfamily_slab,
        size = base_size * 1.7,
        face = "bold",
        lineheight = .8,
        box.color = NA,
        margin = ggplot2::margin(t = 0, b = base_size * .67)
      ),
      plot.subtitle = ggtext::element_textbox_simple(
        family = fontfamily_slab,
        size = base_size,
        lineheight = 1.2,
        color = light_col,
        margin = ggplot2::margin(t = 0, b = base_size * 1.5)
      ),
      plot.caption = ggtext::element_textbox_simple(
        family = fontfamily_slab,
        size = base_size / 2,
        lineheight = 1.2,
        color = light_col,
        margin = ggplot2::margin(t = base_size * 1.5, b = 0)
      ),
      axis.title.x = ggplot2::element_text(
        family = fontfamily_slab,
        margin = ggplot2::margin(t = base_size / 3, r = 3, b = 3, l = 3)
      ),
      axis.title.y = ggplot2::element_text(
        family = fontfamily_slab,
        margin = ggplot2::margin(t = 3, r = 3, b = base_size / 3, l = 3)
      ),
      axis.text.x = ggplot2::element_text(
        color = light_col,
        margin = ggplot2::margin(t = base_size / 4, r = 1, b = 1, l = 1)
      ),
      axis.text.y = ggplot2::element_text(
        color = light_col,
        margin = ggplot2::margin(t = 1, r = base_size / 4, b = 1, l = 1)
      ),
      axis.ticks.length = grid::unit(.33, "lines"),
      strip.text = ggplot2::element_text(
        family = fontfamily_slab,
        face = "bold"
      ),
      strip.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        color = "grey87",
        size = .4
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(
        color = "transparent",
        fill = "transparent"
      ),
      plot.background = ggplot2::element_rect(
        color = "white",
        fill = "white"
      ),
      plot.margin = ggplot2::margin(
        t = margin, r = margin, l = margin, b = margin
      ),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      legend.position = legend_position,
      legend.text = ggplot2::element_text(
        family = fontfamily_slab,
        color = light_col,
        size = base_size * .75
      )
    )

  ## legend title
  if(!legend_title) {
    out <- out + ggplot2::theme(legend.title = ggplot2::element_blank())
  } else {
    out <- out +
      ggplot2::theme(legend.title = ggplot2::element_text(
        family = fontfamily_slab,
        color = light_col,
        size = base_size * .75,
        margin = ggplot2::margin(b = 10)
        )
      )
  }
  
  ## rotate x axis
  if(axis_x_rotate) {
    out <- out +
      ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 1, angle = 45))
  }
  
  ## remove axis titles if selected
  if (axistitle != "none") {
    if (!grepl("X|x", axistitle)) {
      out <- out +
        ggplot2::theme(axis.title.x = ggplot2::element_blank())
    }
    if (!grepl("Y|y", axistitle)) {
      out <- out +
        ggplot2::theme(axis.title.y = ggplot2::element_blank())
    }
  } else {
    out <- out +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank())
  }
  
  ## remove axis text if selected
  if (axistext != "none") {
    if (!grepl("X|x", axistext)) {
      out <- out +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank())
    }
    if (!grepl("Y|y", axistext)) {
      out <- out +
        ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank())
    }
  } else {
    out <- out +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank())
  }
  
  return(out)
}

