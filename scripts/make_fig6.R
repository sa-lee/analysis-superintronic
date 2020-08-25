#' Generate coverage plots for IRfinder and ISA
#' 
#' Generate coverage plots for short and long genes
#' 
#' Note: that the IRfinder results were obtained
#' from github.com/charitylaw/intron-reads
#' 
#' 
#' pretty coverage
library(plyranges)
library(superintronic)
library(ggplot2)
source(here::here("R", "prettycov.R"))

# previously obtained coverage and design
design <- readRDS(here::here("data", "design.rds"))
parts <- readRDS(here::here("data", "complete-annotation.rds"))
cvg <- readRDS(here::here("data", "complete-coverage.rds"))
names(cvg) <- design$Sample

design <- subset(design, Kit == "polyA") 
cvg <- cvg[design$Sample]


targets <- filter(parts, 
                  gene_name %in% c("NBEAL2", "HNRNPL", "HLA-B")) %>% 
  select(gene_id, gene_name) %>% 
  as.data.frame()

# targets <- targets[c(3,1,2),]

.custom_theme <- theme(axis.title.x = element_blank(),
                       axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       strip.text = element_text(hjust = 0, size = 10),
                       axis.text.y = element_text(size = 8),
                       strip.background = element_blank())

38836780- 38837383:-
panel_a <- pretty_cov_plot(cvg, parts, targets[1,], 
                           design, 
                           base_size = 12,
                           cvg_theme = .custom_theme,
                           highlight = as_granges(data.frame(seqnames = "chr19",
                                                             start = 38837586L,
                                                             end = 38837593L)),
                           heights = c(2, 0.25))


panel_b <- pretty_cov_plot(cvg, parts, targets[2,], 
                            design, 
                            base_size = 12,
                            cvg_theme = .custom_theme,
                            highlight = as_granges(data.frame(seqnames = "chr3",
                                                              start = 46996833L,
                                                              end = 46996953L)),
                            heights = c(2, 0.25))

panel_c <- pretty_cov_plot(cvg, parts, targets[3,], 
                           design, 
                           base_size = 12,
                           cvg_theme = .custom_theme,
                           heights = c(2, 0.25))


fig6 <- patchwork::wrap_plots(panel_a, panel_b, panel_c, ncol = 3)

ggsave(here::here("img/Fig6.pdf"), fig6, width = 11)
