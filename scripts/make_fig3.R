#' Coverage over exon/intron regions as histograms
#' Make figure 3 panel a
#' These are line plots, where we have taken absolute 
#' ranges on the x and on the y-axis computed
#' the relative log-coverage (the log-coverage divided
#' by the maximum log-coverage over a gene)
#' These are plotted only on the HCC287 celllines.
library(plyranges)
library(superintronic)
library(ggplot2)

cvg <- readRDS(here::here("data", "parts-coverage.rds")) %>% 
  filter(CellLine == "HCC287") %>%
  mutate(log_score = log2(score + 0.5),
         CellLine = ifelse(CellLine == "HCC287", "HCC827", CellLine)) 

# our selected genes
parts <- readRDS(here::here("data", "filtered-annotation.rds")) 

# filter out low coverage genes
# and add normed score
cvg2 <- cvg %>% 
  group_by(Kit, gene_id) %>% 
  filter(mean(log_score) >= log2(3.5)) %>% 
  ungroup() 


max_cov <- cvg2 %>% 
  group_by(Sample, gene_id) %>% 
  summarise(max_log_score = max(log_score))

# join it back on and add strand information
cvg2<- cvg2 %>% 
  mutate(strand = feature_strand,
         max_log_score = max_cov[match(
           paste0(Sample, gene_id), 
           paste0(max_cov$Sample, max_cov$gene_id)
         ), "max_log_score"],
          normed_score = log_score / max_log_score
  )

# order by coordianates
cvg2 <- sort(cvg2)

# reshape via superintronic
by_gene <- cvg2 %>% 
  group_by(gene_id) %>% 
  dplyr::group_split()


# for each gene section it into twenty bins
library(BiocParallel)
register(MulticoreParam(workers = 4))

bins <- bplapply(by_gene, 
               function(.x) {
                 unlist(
                   GenomicRanges::tile(
                     reduce_ranges_directed(.x), 
                     n = 20
                   )
                 ) %>% 
                   mutate(bin = ifelse(strand == "-", 20:1, 1:20))
               }) %>% 
  as("GRangesList")

# overlap it with genes
olaps <- bplapply(
  seq_along(bins), 
  function(i) join_overlap_intersect_directed(bins[[i]], by_gene[[i]])
  ) %>% 
  as("GRangesList")


olaps <- unlist(olaps) %>% 
  mutate(bin = Rle(bin))


# summarise over bins for each gene and feature
bin_means <- olaps %>% 
  group_by(Sample, gene_id, feature_type, bin) %>% 
  summarise(mn = mean(normed_score)) %>% 
  dplyr::as_tibble() 

gene_lengths <- parts %>% 
  select(gene_id, width, .drop_ranges = TRUE) %>% 
  dplyr::as_tibble() %>% 
  mutate(
    cat = dplyr::case_when(
      width >= quantile(width, 2/3) & width <= max(width) ~ "long",
      width >= quantile(width, 1/3) & width < quantile(width, 2/3) ~ "regular",
      width >= min(width) & width < quantile(width, 1/3) ~ "short"
    ),
    cat = factor(cat, levels = c("short", "regular", "long"))
  )

# summarise over all bins within each sample and feature
tbl <- bin_means %>% 
  dplyr::left_join(gene_lengths) %>% 
  group_by(Sample, cat, feature_type, bin) %>% 
  summarise(lumpy = mean(mn), bumpy = var(mn)) %>% 
  ungroup() %>% 
  mutate(Kit = gsub("R[1-3]_HCC287_", "", Sample),
         feature_type = factor(feature_type, levels = c("exon", "intron"))) 

fig3a <- ggplot(tbl, aes(x = bin, y = lumpy, colour = Kit, group = Sample)) +
  geom_line(size = 1) +
  facet_grid(feature_type ~ cat, scales = "free_y") +
  scale_x_continuous(limits = c(1,20)) +
  scale_color_manual(values = c("#404040", "#bababa")) +
  theme_bw(base_size = 20) +
  theme(aspect.ratio = 1,
        legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid = element_blank()
        ) +
  labs(y = "Relative Log-coverage", 
       x = "Bin (Oriented 5' to 3')")

ggsave(here::here("img", "Fig3a.pdf"), fig3a, width = 11)  


source(here::here("R", "prettycov.R"))

# replicate figure 3b, two short two long
gene_names <- c("ENSG00000136997.17", "ENSG00000117525.13", 
                "ENSG00000196937.10", "ENSG00000196428.12")
names(gene_names) <- paste0("(", letters[1:4], ")")

targets <- filter(parts, 
                  gene_id %in% c("ENSG00000136997.17", "ENSG00000117525.13", 
                                 "ENSG00000196937.10", "ENSG00000196428.12")
)

.custom_theme <- theme(axis.title.x = element_blank(),
                       axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       strip.text = element_text(hjust = 0, size = 10),
                       axis.text.y = element_text(size = 8),
                       strip.background = element_blank())

short <- lapply(1:2, 
                function(.x) pretty_cov_plot(cvg, parts, 
                                             targets[.x,], 
                                             alpha = TRUE, 
                                             base_size = 12,
                                             cvg_theme = .custom_theme,
                                             .label = names(gene_names)[.x],
                                             heights = c(2, 0.25))) %>% 
  patchwork::wrap_plots(ncol = 2, guides = "keep")

long <- lapply(3:4, function(.x) pretty_cov_plot(cvg, parts, 
                                                 targets[.x,], 
                                                 alpha = TRUE, 
                                                 base_size = 12,
                                                 cvg_theme = .custom_theme,
                                                 .label = names(gene_names)[.x],
                                                 heights = c(2, 0.25))) %>% 
  patchwork::wrap_plots(ncol = 2, guides = "keep")

fig3b <- patchwork::wrap_plots(short, long, nrow = 2, guides = "keep")

ggsave(here::here("img", "Fig3b.pdf"), fig3b, width = 11)  

# fig3 <- patchwork::wrap_plots(fig3a, fig3b, ncol = 2, guides = "keep",
#                               widths = c(3,4.75))
# 
# ggsave(here::here("img", "Fig3.png"), fig3, width = 16, height = 8, 
#        dpi = "retina")
