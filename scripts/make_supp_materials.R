#' Make supplementary figure of all 
#' - superintronic hits 
#' - ir finder hits

#' superintronic hits

library(plyranges)
library(superintronic)
library(ggplot2)
source(here::here("R", "prettycov.R"))

cvg <- readRDS(here::here("data", "parts-coverage.rds"))

parts <- readRDS(here::here("data", "complete-annotation.rds"))

super <- readRDS(here::here("data", "superintronic-hits.rds")) %>% 
  filter(CellLine == "HCC287", Kit == "polyA") %>% 
  ungroup()

cvg_sub <- filter(cvg, Kit == "polyA") %>% 
  mutate(CellLine = ifelse(CellLine == "HCC287", "HCC827", CellLine)) %>% 
  mutate(log_score = log2(score + 0.5))

for (i in seq_len(nrow(super))) {
  input <- super[i,]
  track_plot <- pretty_cov_plot(cvg_sub, parts, input, heights = c(2,0.25))
  fn <- here::here("img", "superintronic-polyA-cov", 
                   paste0("superintronic-polyA-", input$gene_id, ".pdf"))
  ggsave(fn, track_plot)
}

#' ir finder hits
library(dplyr)
design <- readRDS(here::here("data", "design.rds"))
cvg <- readRDS(here::here("data", "complete-coverage.rds"))
names(cvg) <- design$Sample

design <- subset(design,  Kit == "polyA") 

cvg_sub <- cvg[design$Sample]


ir_finder_ac <- readr::read_tsv(here::here("data", "IRfinder_AC_hits.tsv"), 
                                skip = 8,
                                col_types = c("Chr" = "c", 
                                              "Start" = "i", 
                                              "End" = "i")) %>% 
  rename_all(.funs = ~stringr::str_replace_all(stringr::str_to_lower(.), "-|\\/", "_")) %>% 
  filter(p_diff < 33) %>% 
  mutate(splitter = stringr::str_split(intron_genename_geneid,  "\\/"),
         gene_name = purrr::map_chr(splitter, ~.[[1]]),
         gene_id = purrr::map_chr(splitter, ~.[[2]])
  ) %>% 
  select(chr, start, end, gene_id, gene_name, p_diff)


parts_tbl <- parts %>% 
  plyranges::select(gene_id, gene_name, .drop_ranges = TRUE) %>% 
  as_tibble()

irf_ac <- ir_finder_ac %>% 
  left_join(parts_tbl, by = "gene_name") %>% 
  select(gene_id = gene_id.y, gene_name, adj_pval = p_diff) %>% 
  mutate(`IRFinder_AC` = 1L) %>% 
  as.data.frame()

for (i in seq_len(nrow(irf_ac))) {
  input <- irf_ac[i,]
  track_plot <- pretty_cov_plot(cvg_sub, parts, input, design, heights = c(2,0.25))
  fn <- here::here("img", "irfinder-cov", 
                   paste0("ac-test", input$gene_id, ".pdf"))
  ggsave(fn, track_plot)
}


ir_finder_glm <- readr::read_tsv(here::here("data", "IRfinder_GLM_hits.tsv"),
                                 skip = 1,
                                 col_names = c("id",
                                               "base_mean",
                                               "logFC",
                                               "logFC_se",
                                               "stat",
                                               "pval",
                                               "adj_pval")) %>% 
  filter(adj_pval < 0.05) %>% 
  arrange(adj_pval) %>% 
  mutate(splitter = stringr::str_split(id,  "\\/"),
         gene_name = purrr::map_chr(splitter, ~.[[1]]),
         gene_id = purrr::map_chr(splitter, ~.[[2]]),
         region = stringr::str_replace_all(purrr::map_chr(splitter, ~.[[4]]), " ", ""),
  ) %>% 
  select(gene_name, gene_id, region, adj_pval)

irf_glm <- ir_finder_glm %>% 
  group_by(gene_name, gene_id) %>% 
  tidyr::nest() %>% 
  left_join(parts_tbl, by = "gene_name") %>% 
  left_join(parts_tbl %>% 
              mutate(gene_id2 = stringr::str_replace(gene_id, "\\.[1-9]{1,}", "")), 
            by = c("gene_id.x" = "gene_id2")) %>% 
  group_by(gene_name.x) %>% 
  mutate(canon_gene_id = list(unique(c(gene_id, gene_id.y))),
         canon_gene_id = lapply(canon_gene_id, function(x) {
           if (any(!is.na(x))) return(x[!is.na(x)])
           return(NA_character_)
         }),
         canon_gene_id = as.character(unlist(canon_gene_id))) %>% 
  ungroup() %>% 
  select(gene_name = gene_name.y, gene_id = canon_gene_id, data) %>% 
  filter(!is.na(gene_id)) %>% 
  mutate(`IRFinder_GLM` = 1L)

for (i in seq_len(20)) {
  input <- irf_glm[i,]
  track_plot <- pretty_cov_plot(cvg_sub, parts, input, design, heights = c(2,0.25))
  fn <- here::here("img", "irfinder-cov", 
                   paste0("glm-test", input$gene_id, ".pdf"))
  ggsave(fn, track_plot)
}
