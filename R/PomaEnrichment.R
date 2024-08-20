
#' Enrichment Analysis
#'
#' @description `PomaEnrichment` performs enrichment analysis on a set of query gene symbols using specified methods and gene set collections. It allows for the analysis of over-representation (ORA) or gene set enrichment (GSEA) in various model organisms.
#'
#' @param genes Character vector. Set of query gene symbols.
#' @param method Character. Enrichment method. Options are: 'ora' (simple over-representation analysis based on hypergeometric test) and 'gsea' (gene set enrichment analysis on a ranked list of genes).
#' @param organism Character. Indicates the model organism name. Default is 'Homo sapiens'. Other options are: 'Anolis carolinensis', 'Bos taurus', 'Caenorhabditis elegans', 'Canis lupus familiaris', 'Danio rerio', 'Drosophila melanogaster', 'Equus caballus', 'Felis catus', 'Gallus gallus', 'Macaca mulatta', 'Monodelphis domestica', 'Mus musculus', 'Ornithorhynchus anatinus', 'Pan troglodytes', 'Rattus norvegicus', 'Saccharomyces cerevisiae', 'Schizosaccharomyces pombe 972h-', 'Sus scrofa', 'Xenopus tropicalis'. See `msigdbr::msigdbr_show_species()`.
#' @param collection Character. Indicates the gene set collection. Default is 'C5' (Gene Ontology gene sets). Other options are: 'C1' (positional gene sets), 'C2' (curated gene sets), 'C3' (regulatory target gene sets), 'C4' (computational gene sets), 'C6' (oncogenic signature gene sets), 'C7' (immunologic signature gene sets), 'C8' (cell type signature gene sets), 'H' (Hallmark gene sets). See `msigdbr::msigdbr_collections()`.
#' @param universe Character vector. A universe from which 'genes' were selected.
#' @param rank Numeric vector. Ranking factor to sort genes for GSEA (e.g., logFC, -log10(p-value), etc).
#' @param pval_cutoff Numeric. Raw p-value cutoff on enrichment tests to report.
#' @param fdr_cutoff Numeric. Adjusted p-value cutoff on enrichment tests to report.
#' @param min_size Numeric. Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param max_size Numeric. Maximal size of a gene set to test. All pathways above the threshold are excluded.
#'
#' @export
#'
#' @return A `tibble` with the enriched gene sets.
#' @author Pol Castellano-Escuder
#'
#' @importFrom magrittr %>%
#' 
#' @examples
#' # Example genes
#' genes <- c("BRCA1", "TP53", "EGFR", "MYC", "PTEN")
#' 
#' # Perform ORA on Gene Ontology (C5) gene sets for Homo sapiens
#' PomaEnrichment(
#'   genes = genes,
#'   method = "ora",
#'   organism = "Homo sapiens",
#'   collection = "C5",
#'   pval_cutoff = 0.05,
#'   fdr_cutoff = 0.1,
#'   min_size = 10,
#'   max_size = 500)
#' 
#' # Example genes with ranking factors (e.g., logFC values)
#' genes <- c("Actb", "Gapdh", "Cdkn1a", "Cd44", "Pten")
#' rank <- c(2.5, -1.8, 3.1, -2.2, 1.7)
#' 
#' # Perform GSEA on Hallmark (H) gene sets for Mus musculus
#' PomaEnrichment(
#'   genes = genes,
#'   method = "gsea",
#'   organism = "Mus musculus",
#'   collection = "H",
#'   rank = rank,
#'   pval_cutoff = 0.05,
#'   fdr_cutoff = 0.25,
#'   min_size = 15,
#'   max_size = 500)
PomaEnrichment <- function(genes,
                           method = "ora",
                           organism = "Homo sapiens", # msigdbr::msigdbr_show_species()
                           collection = "C5", # msigdbr::msigdbr_collections()
                           universe = NULL,
                           rank = NULL,
                           pval_cutoff = 0.05,
                           fdr_cutoff = 0.1,
                           min_size = 2,
                           max_size = if (method == "gsea") {length(genes) - 1} else {NULL}) {
  
  pathways <- msigdbr::msigdbr(species = organism, 
                               category = collection,
                               subcategory = NULL)
  
  pathways_format <- split(x = pathways$gene_symbol, f = pathways$gs_name)

  if (method == "gsea") {
    if (is.null(rank)) {
      stop("ranking variable is missing")
    }
    
    names(rank) <- genes
    ordered_genes <- sort(rank, decreasing = TRUE)
    
    gsea_result <- fgsea::fgsea(pathways = pathways_format,
                                stats = ordered_genes,
                                eps = 0,
                                minSize = min_size,
                                maxSize = max_size,
                                scoreType = "std")
    
    result <-
      gsea_result %>%
      dplyr::rowwise() %>%
      dplyr::mutate(leadingEdge = paste0(leadingEdge, collapse = ", ")) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(direction = ifelse(NES > 0, "Up", "Down")) %>%
      dplyr::select(pathway, NES, pval, adjPval = padj,
                    direction, pathway_size = size, leading_edge = leadingEdge) %>%
      dplyr::arrange(pval) %>%
      dplyr::as_tibble() %>%
      dplyr::filter(pval < pval_cutoff & adjPval < fdr_cutoff)
  }
  
  else if (method == "ora") {
    if (is.null(universe)) {
      universe <- unique(pathways$gene_symbol)
    }
    if (is.null(max_size)) {
      max_size <- length(unique(pathways$gene_symbol)) - 1
    }
    
    suppressWarnings({
      ora_result <- fgsea::fora(pathways = pathways_format,
                                genes = genes,
                                universe = universe,
                                minSize = min_size,
                                maxSize = max_size)
    })
    
    result <-
      ora_result %>%
      dplyr::rowwise() %>%
      dplyr::mutate(overlapGenes = paste0(overlapGenes, collapse = ", ")) %>%
      dplyr::ungroup() %>%
      dplyr::select(pathway, pval, adjPval = padj,
                    pathway_size = size, overlap, overlap_genes = overlapGenes) %>%
      dplyr::arrange(pval) %>%
      dplyr::as_tibble() %>%
      dplyr::filter(pval < pval_cutoff & adjPval < fdr_cutoff)
  }
  
  return(result)
}

