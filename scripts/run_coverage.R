# Batch Process - coverage analysis
# This runs the 'superintronic' pipeline for 
# computing coverage over the intronic/exonic parts of a gene
library(here)
library(BiocParallel)
suppressPackageStartupMessages(library(plyranges))
suppressPackageStartupMessages(library(superintronic))

# --- design table ---
design <- read.table(
  here("data-raw", "targets.txt"),
  header = TRUE,
  stringsAsFactors = FALSE
) %>% 
  S4Vectors::DataFrame() %>% 
  transform(
    File = S4Vectors::Rle(here("data-raw", "bam", sub("\\.", "-", File)))
  ) %>% 
  BiocGenerics::subset(Replicate %in% c("R1", "R2", "R3") & 
                         Mixture %in% c(0,100),
                       select = c(File, Mixture, Replicate, Kit)) %>%  
  S4Vectors::transform(
    CellLine = ifelse(Mixture == 0, "HCC287","NCI-H1975"),
    Kit = S4Vectors::Rle(ifelse(Kit == "mRNA", "polyA", "total-RNA"))
  ) %>% 
  S4Vectors::transform(
    Sample =  paste(Replicate, CellLine, Kit, sep = "_")
)

print(design)

saveRDS(design, here("data", "design.rds"))


# --- annotation ---
gff <- here("data-raw", "gencode.v27.annotation.gtf.gz")

gr_gff <- read_gff(gff, genome_info = "hg38")
parts <- collect_parts(gr_gff)

print(head(parts))

saveRDS(parts, here("data", "complete-annotation.rds"))

# --- filter annotation ---
parts_sub <- parts %>% 
  filter(
    gene_type == "protein_coding", 
    n_olaps == 1, 
    seqnames != "chrM",
    lengths(exonic_parts) > 1
  ) %>%
  GenomeInfoDb::keepStandardChromosomes("Homo sapiens", "coarse") %>% 
  GenomeInfoDb::dropSeqlevels("chrM", "coarse")

print(parts_sub)


saveRDS(parts_sub, here("data", "filtered-annotation.rds"))

# --- compute coverage ---

# store it as a lits for space purposes 
bams <- Rsamtools::BamFileList(as.character(design$File))

info <-  bams %>% 
  get_genome_info() %>% 
  GenomeInfoDb::keepStandardChromosomes("Homo sapiens", "coarse") %>% 
  set_genome_info(genome = "hg38")


cvg_rng <- bplapply(bams, 
                    function(.) filter_by_overlaps(compute_coverage(.), info),
                    BPPARAM = MulticoreParam(4))

cvg_rng <- GRangesList(cvg_rng)

names(cvg_rng) <- design$Sample 

saveRDS(cvg_rng, here("data", "complete-coverage.rds"))


# --- merge parts ---
cvg_over_features <- bplapply(cvg_rng,
                              function(x) join_parts(x, parts_sub),
                              BPPARAM =  MulticoreParam(4))
cvg_over_features <- as(cvg_over_features, "GRangesList")

md <- DataFrame(Sample = Rle(names(cvg_over_features), lengths(cvg_over_features)),
                CellLine = Rle(design$CellLine, lengths(cvg_over_features)),
                Kit = Rle(as.factor(design$Kit), lengths(cvg_over_features)),
                Replicate = Rle(design$Replicate, lengths(cvg_over_features)))

cvg_gr <- unlist(cvg_over_features, use.names = FALSE)

mcols(cvg_gr) <- cbind(mcols(cvg_gr), md)

saveRDS(cvg_gr, here("data", "parts-coverage.rds"))
