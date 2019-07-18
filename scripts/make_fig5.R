#' Figure 5:
#' Generate coverage plots for three genes that 
#' were superintronic hits: ARGLU1, PSMB7 EIF2S3
#' These are plotted only on the HCC287 celllines.
library(superintronic)
library(plyranges)
library(ggplot2)
library(patchwork)

coverage_plot <- function(target) {
  cvg_rng <- filter_by_overlaps(cvg, target) %>% 
    mutate(strand = strand(target),
           var = paste0("Kit: ", Kit, " Cellline: ", CellLine))
  p <- view_coverage(cvg_rng, 
                     score = log_score, 
                     colour = feature_type,  
                     facets = dplyr::vars(var)) +
    ylim(0, NA) +
    guides(colour = FALSE) +
    scale_color_brewer(palette = "Dark2") +
    labs(subtitle = paste("Coverage over", target$gene_name),
         y = "Log-coverage") +
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

parts_sub <- readRDS(here::here("data", "filtered-annotation.rds"))

cvg <- readRDS(here::here("data", "parts-coverage.rds")) %>% 
  filter(CellLine == "HCC287") %>% 
  mutate(log_score = log2(score + 0.5))


gene_names <- c("PSMB7", "EIF2S3", "ARGLU1")


plot_lists <- lapply(gene_names, function(.x) {
  target <- filter(parts_sub, gene_name == !!.x)
  coverage_plot(target)
})

fig <- wrap_plots(plot_lists, ncol = 1, guides = "keep")

ggsave(here::here("img/Fig5.pdf"), fig, height = 7, width = 12)
