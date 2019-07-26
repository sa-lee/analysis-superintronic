#' Pretty track plot
#' 
#' @param cvg Coverage as GRanges or GRangesList
#' @param parts Annotation from `superintronic::collect_parts()`
#' @param target Region of interest as GRanges or data.frame
#' @param alpha (Make low coverage regions transparent?) Default is FALSE
#' @param ... Other options passed to `patchwork::wrap_plots()`
#' 
#' @import ggplot2 GenomicRanges plyranges superintronic S4Vectors
#' @importFrom dplyr case_when
#' @importFrom methods is
#' @importFrom patchwork wrap_plots
#' @importFrom BiocGenerics strand mean unique
#' @export
pretty_cov_plot <- function(cvg, parts, target, design = NULL, alpha = FALSE, ...) {
  
  if(is(target, "data.frame")) {
    target <- plyranges::filter(parts, gene_id == !!target$gene_id)
  }
  
  cvg_rng <- set_plot_rng(cvg, target, design, alpha)
  
  if (alpha) {
    p <- ggplot(as.data.frame(cvg_rng), 
                aes(x = start, xend = end, y = 0, yend = log_score)) +
      geom_segment(aes(colour = feature_type, alpha = alpha)) +
      superintronic:::rescale_by_width(cvg_rng) +
      guides(alpha = FALSE)
  } else {
    p <- superintronic::view_coverage(cvg_rng, 
                                      score = log_score, 
                                      colour = feature_type,  
                                      facets = dplyr::vars(var))
  }
  
  p <- p +     
    ylim(0, NA) +
    guides(colour = FALSE) +
    scale_color_brewer(palette = "Dark2") +
    labs(subtitle = paste("Coverage over", target$gene_name),
         y = "Log coverage"
    ) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text = element_text(hjust = 0, size = 8),
          strip.background = element_blank()) +
    expand_limits(y = 0)
  
  segments <- superintronic::unnest_parts(target)
  track <- view_segments(segments, colour = feature_type) +
    theme(axis.text.x = element_text(size = 8))
  
  patchwork::wrap_plots(p, 
                        track, 
                        ncol = 1, 
                        ...)
  
}



set_plot_rng <- function(cvg, target, design, alpha) {
  if (is(cvg, "GRangesList")) {
    stopifnot(is(design, "DataFrame"))
    cvg_rng <- as(lapply(cvg, function(x) superintronic::join_parts(x, target)), 
                  "GRangesList")
    md <- S4Vectors::DataFrame(
      Sample = S4Vectors::Rle(names(cvg_rng), lengths(cvg_rng)),
      CellLine = S4Vectors::Rle(design$CellLine, lengths(cvg_rng)),
      Kit = S4Vectors::Rle(as.factor(design$Kit), lengths(cvg_rng))
    )
    
    cvg_rng <- unlist(cvg_rng, use.names = FALSE)
    S4Vectors::mcols(cvg_rng) <- cbind(S4Vectors::mcols(cvg_rng), md)
    
    cvg_rng <- plyranges::mutate(cvg_rng,
                                 log_score = log2(score + 0.5),
                                 strand = BiocGenerics::strand(target),
                                 var = paste0("Kit: ", Kit, " Cellline: ", CellLine))
  }
  
  if (is(cvg, "GRanges")) {
    cvg_rng <- plyranges::filter_by_overlaps(cvg, target)
    cvg_rng <- plyranges::mutate(cvg_rng, 
                                 strand = BiocGenerics::strand(target),
                                 var = paste0("Kit: ", Kit, " Cellline: ", CellLine))
  }
  
  if (alpha) {
    cvg_rng <- plyranges::disjoin_ranges_directed(group_by(cvg_rng, var),
                                                  log_score = BiocGenerics::mean(log_score),
                                                  feature_type = unlist(BiocGenerics::unique(feature_type)))
    
    cvg_rng <- plyranges::mutate(cvg_rng,
                                 alpha = dplyr::case_when(
                                   score > 3 & feature == "intron" ~ 1,
                                   score <= 3 & feature == "intron"~ 0.5,
                                   TRUE ~ 1
                                 ))
  }
  
  return(cvg_rng)
}