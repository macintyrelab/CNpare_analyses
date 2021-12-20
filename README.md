# CNpare_analyses

## Summary of the repository 
This repository includes the code needed to reproduce the analyses performed in the manuscript [CNpare: matching DNA copy number profiles](https://www.biorxiv.org/content/10.1101/2021.09.28.462193v1) 

Prior to run this code, our **CNpare** tool should be installed as follows:
``` r
git clone https://github.com/macintyrelab/CNpare.git
```

Once installed, CNpare must be loaded
``` r
library(CNpare)
```

## Source code

The source code for reproducing all analysess is split across four files:

* `CNpare_AssessingTool.Rmd` -- This Rmarkdown contains the code for assessing the performance of CNpare by using separate cultures of the same cell lines profiled as part of the CCLE and GDSC projects (304 pairs)
* `Comparing_OtherData.Rmd` -- This Rmarkdown shows the benchmarking of CNpare against other approaches for choosing matching cell lines. These approaches are: gene-level copy number, chromosome arm-level copy number, ploidy status and gene expression. 
* `CNpare_DemonstratingUtility.Rmd` -- This Rmarkdown is an example of the CNpare workflow. We use the high-grade serous ovarian cancer cell line OVKATE to find the next most representative cell line model in the CNpare database. Matching is perfomed by similarity in both copy-number profiles and exposure to copy number signatures
* `Suitability_OVKATEMatches.Rmd` -- This Rmarkdown includes the validation analysis of the suitability of the cell lines matched with OVKATE. To this, we assessed correlation of gene expression between the cell lines matched. 

Data needed for reproducing the analysis is included in this repository, or included in the CNpare package, or can be downloaded from public databases (links included in the Rmarkdowns)
