# Intron retention analysis of human celline mixtures

A compendium of scripts and analysis for Lee, Zhang et. al, 2019 
using `superintronic`, including supplementary information.

Other scripts and analyses are available as part of https://github.com/charitylaw/intron-reads

The project is laid out as an R package using the following directory structure:



- `data/` (cached data, saved from Rmd and R scripts, includes kallisto results)
- `Rmd/` (reports for each part of the analysis)
- `R/` (functions used to create "pretty" coverage plots)
- `img/` (all figures created during analyses; including supp figures)
- `scripts/` (batch processing of raw data files, and scripts to reproduce figures)


You can install all the analysis dependencies using:

```{r}
# install.packages("BiocManager")
BiocManager::install("sa-lee/analysis-superintronic")
```


Note that you will need to have [git-lfs](https://git-lfs.github.com/)
installed if you would like to use the cached data directly. 

