#' Generate coverage plots for IRfinder and ISA
#' 
#' Note: that the IRfinder results were obtained
#' from github.com/charitlaw/intron-reads
#' 
#' pretty coverage
coverage_plot <- function(.y) {
  target <- filter(parts, gene_id == !!.y$gene_id)
  
  # set up parts
  cvg_rng <- as(lapply(cvg, function(x) join_parts(x, target)), "GRangesList")
  md <- DataFrame(Sample = Rle(names(cvg_rng), lengths(cvg_rng)),
                  CellLine = Rle(design$CellLine, lengths(cvg_rng)),
                  Kit = Rle(as.factor(design$Kit), lengths(cvg_rng)))
  cvg_rng <- unlist(cvg_rng, use.names = FALSE)
  mcols(cvg_rng) <- cbind(mcols(cvg_rng), md)
  cvg_rng <- cvg_rng %>% 
    mutate(log_score = log2(score + 0.5),
           strand = strand(target),
           var = paste0("Kit: ", Kit, " Cellline: ", CellLine))
  
  
  p <- view_coverage(cvg_rng, 
                     score = log_score, 
                     colour = feature_type,  
                     facets = dplyr::vars(var)) +
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
  
  segments <- superintronic:::flatten_parts(target)
  track <- view_segments(segments, color = feature_type) +
    theme(axis.text.x = element_text(size = 8))
  
  patchwork::wrap_plots(p, track, 
                        ncol = 1, 
                        heights = c(2, 0.25))
  
}


library(plyranges)
library(superintronic)


# previously obtained coverage and design
design <- readRDS(here::here("data", "design.rds"))
parts <- readRDS(here::here("data", "complete-annotation.rds"))
cvg <- readRDS(here::here("data", "complete-coverage.rds"))

design <- subset(design, Kit == "polyA")
cvg <- cvg[design$Sample]


# load in IRFinder results



targets <- filter(parts, 
                  gene_name %in% c("NBEAL2", "HNRNPL", "HLA-B")
) %>% 
  select(gene_id) %>% 
  as.data.frame()

fig6_plots <- lapply(seq_len(nrow(targets)),
                     function(.x) coverage_plot(targets[.x,]))

fig6 <- patchwork::wrap_plots(fig6_plots, ncol = 1, guides = "keep")

ggsave(here::here("figures/Fig6.pdf"), fig6, height = 7, width = 12)

