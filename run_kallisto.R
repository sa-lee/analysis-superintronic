#!/usr/bin/env Rscript
# run kallisto/salmon on targets
library(dplyr, warn.conflicts = FALSE)

find_files <- function(x, path) {
  files <- list.files(path, full.names = TRUE)
  files[grep(x, files)]
}

# fastq path
fastq_path <- file.path("fastq")
# index
index <- here::here("data", "transcripts.idx")

# read in design table
samples <- read.table(
  here::here("data-raw", "targets.txt"),
  header = TRUE,
  stringsAsFactors = FALSE
  ) %>% 
  as_tibble() %>% 
  filter(Mixture == 0, Replicate %in% c("R1", "R2", "R3")) %>% 
  mutate(fq_target = sub("_R1.bam", "", basename(File))) %>%
  mutate(fq_target = sub("\\.", "-", fq_target)) %>% 
  mutate(output = paste0(sub("\\.", "-", Sample), "_", Kit)) %>% 
  mutate(fq_files = lapply(fq_target, function(x) find_files(x, fastq_path))) %>% 
  select(output, fq_target, fq_files)
  
# ready to run!
# setup commands

system("module load kallisto/0.44.0")


command <- "kallisto quant"
args <- c("-i", index, "-o",  NA, 
          "--single", "-l", 200, "-s", 20)

# loop over each sample
for (i in seq_along(samples)) {
  row <- samples[i, ]
  cat(paste("Preparing salmon call for", row$output), "\n")
  outdir <- here::here("data", "kallisto", row$output)
  
  if (!dir.exists(outdir)) dir.create(outdir)
  
  args[4] <- outdir
  
  fastqs <- paste(row$fq_files[[1]], collapse = " ")
  call <- paste(command, paste(args, collapse = " "), fastqs)
  cat(call, "\n")
  system(call)
  
}