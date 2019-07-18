# Intron retention analysis of human celline mixtures

A compendium of scripts and analysis for Lee, Zhang et. al, 2019 (preprint link coming) using `superintronic`, and `IsoformSwitchAnalyzeR`. Other scripts and analyses 
are available as part of https://github.com/charitylaw/intron-reads

The project is laid out as an R package using the following directory structure:

```{r}
data-raw/ (raw data for celllines: fastq files, bam files and reference; not included with repository)
data/ (cached data, saved from Rmd and R scripts, includes kallisto results)
Rmd/ (reports for each part of the analysis)
img/ (saved figures; including supp figures)
scripts/ (batch processing of data-raw, plus scripts to reproduce figures)
.*ignore (usual git stuff)
.here (tracking for relative paths using the `here` pacakge)
DESCRIPTION
```


You can install all the analysis package 
dependencies using:

```{r}
# install.packages("BiocManager")
BiocManager::install("sa-lee/analysis-superintronic")
```


Note that you will need to [git-lfs](https://git-lfs.github.com/)
installed if you would like to use the cached data directly. 

