
#' Volcano Plot
#'
#' @description PomaVolcano() generates a volcano plot from the PomaUnivariate(method = "ttest") result. The data can't have negative values.
#'
#' @param data A SummarizedExperiment object.
#' @param pval Select a pvalue type to generate the volcano plot. Options are: "raw" and "adjusted".
#' @param pval_cutoff Numeric. Define the pvalue cutoff (horizontal line).
#' @param adjust Multiple comparisons correction method for t test result. Options are: "fdr", "holm", "hochberg", "hommel", "bonferroni", "BH" and "BY".
#' @param log2FC Numeric. Define the log2 fold change cutoff (vertical lines).
#' @param xlim Numeric. Define the limits for x axis.
#' @param labels Logical that indicates if selected labels will be plotted or not. Defaul is FALSE.
#' @param paired Logical that indicates if the data is paired or not.
#' @param var_equal Logical that indicates if the data variance is equal or not.
#' @param interactive Logical that indicates if an interactive plot will be plotted or not. Defaul is FALSE.
#' @param plot_title Logical that indicates if title will be plotted or not. Defaul is TRUE.
#' 
#' @export
#'
#' @return A ggplot2 object.
#' @author Pol Castellano-Escuder
#'
#' @importFrom magrittr %>% %<>%
#' 
#' @examples 
#' data("st000336")
#' 
#' st000336 %>% 
#'   PomaImpute() %>%
#'   PomaVolcano()
PomaVolcano <- function(data,
                        pval = "raw",
                        pval_cutoff = 0.05,
                        adjust = "fdr",
                        log2FC = 0.6,
                        xlim = 2,
                        labels = FALSE,
                        paired = FALSE,
                        var_equal = FALSE,
                        interactive = FALSE,
                        plot_title = TRUE){

  if (missing(data)) {
    stop("data argument is empty!")
  }
  if(!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaSummarizedExperiment or SummarizedExperiment::SummarizedExperiment")
  }
  if (length(table(SummarizedExperiment::colData(data)[,1])) > 2) {
    stop("Your data have more than two groups!")
  }
  if (!(pval %in% c("raw", "adjusted"))) {
    stop("Incorrect value for pval argument!")
  }
  if (!(adjust %in% c("fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY"))) {
    stop("Incorrect value for adjust argument!")
  }

  df <- POMA::PomaUnivariate(data, method = "ttest", adjust = adjust, paired = paired, var_equal = var_equal)

  names <- rownames(SummarizedExperiment::assay(data))

  if(pval == "raw"){
    df <- data.frame(pvalue = df$pvalue, FC = log2(df$FC), names = names)
  }
  else {
    df <- data.frame(pvalue = df$pvalueAdj, FC = log2(df$FC), names = names)
  }

  df %<>%
    dplyr::mutate(threshold = dplyr::case_when(df$pvalue >= pval_cutoff ~ "None",
                                               df$FC < -log2FC ~ "Down-regulated",
                                               df$FC > log2FC ~ "Up-regulated"
                                               )
    )
  
  volcanoP <- ggplot2::ggplot(data = df, ggplot2::aes(x = FC, y = -log10(pvalue), colour = threshold, label = names)) +
    ggplot2::geom_point(size = 1.75) +
    ggplot2::xlim(c(-(xlim), xlim)) +
    ggplot2::labs(x = "log2 Fold Change",
                  y = "-log10 p-value",
                  color = NULL) +
    ggplot2::scale_y_continuous(trans = "log1p")+
    {if(plot_title)ggplot2::ggtitle(paste0(names(table(SummarizedExperiment::colData(data)[,1]))[2], "/", 
                                           names(table(SummarizedExperiment::colData(data)[,1]))[1]))} +
    ggplot2::geom_vline(xintercept = -log2FC, colour = "black", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = log2FC, colour = "black", linetype = "dashed") +
    ggplot2::geom_hline(yintercept = -log10(pval_cutoff), colour = "black", linetype = "dashed") +
    {if(labels)ggrepel::geom_label_repel(data = df[df$pvalue < pval_cutoff & (df$FC > log2FC | df$FC < -log2FC),],
                                         ggplot2::aes(x = FC, y = -log10(pvalue), label = names), show.legend = FALSE)} +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "right") +
    ggplot2::scale_color_manual(values = c("Down-regulated" = "#E64B35", 
                                           "Up-regulated" = "#3182bd", 
                                           "None" = "#636363"))

  if(interactive) {
    if(!(require("plotly", character.only = TRUE))){
      warning("Package 'plotly' is required for an interactive volcano plot\nUse 'install.packages('plotly')'")
    } else {
      volcanoP <- plotly::ggplotly(volcanoP)
      }
    }

  return(volcanoP)
  
}

