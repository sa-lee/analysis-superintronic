#' Figure 5:
#' Generate coverage plots for three genes that 
#' were superintronic hits: ARGLU1, PSMB7 EIF2S3
#' These are plotted only on the HCC287 celllines.
library(superintronic)
library(plyranges)
library(ggplot2)
library(patchwork)

source(here::here("R", "prettycov.R"))

parts_sub <- readRDS(here::here("data", "filtered-annotation.rds"))

cvg <- readRDS(here::here("data", "parts-coverage.rds")) %>% 
  filter(CellLine == "HCC287") %>% 
  mutate(log_score = log2(score + 0.5),
         CellLine = ifelse(CellLine == "HCC287", "HCC827", CellLine))


gene_names <- c("PSMB7", "EIF2S3", "ARGLU1")


plot_lists <- lapply(gene_names, function(.x) {
  target <- filter(parts_sub, gene_name == !!.x)
  pretty_cov_plot(cvg, parts_sub, target, heights = c(2, 0.25))
})

fig <- wrap_plots(plot_lists, ncol = 1, guides = "keep")

ggsave(here::here("img/Fig5.pdf"), fig, height = 6, width = 10)
