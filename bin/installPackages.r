setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")
set.seed(11)

cat(":: Installing / loading the R packages ::\n:: required by the script ::\n\n")

# Here we install the CRAN missing packages
#  "markdown", "knitr", "rmarkdown", "pacman", 
list.of.packages <- c("easypackages", "geneExpressionFromGEO","dplyr", "jetset") # other packages
new_packages_to_install <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new_packages_to_install)) install.packages(new_packages_to_install, repos="https://utstat.toronto.edu/cran/")

# Here we install the Bioconductor missing packages

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos="https://utstat.toronto.edu/cran/")
#    "annotate"
listOfBiocPackages <- c("Biobase", "org.Hs.eg.db",  "GEOquery", "diffcoexp", "annotate", "sva")

bioCpackagesNotInstalled <- which( !listOfBiocPackages %in% rownames(installed.packages()) )
if(length(bioCpackagesNotInstalled)) cat("package missing listOfBiocPackages[", bioCpackagesNotInstalled, "]: ", listOfBiocPackages[bioCpackagesNotInstalled], "\n", sep="")

# check there's still something left to install
if( length(bioCpackagesNotInstalled) ) {
    BiocManager::install(listOfBiocPackages[bioCpackagesNotInstalled])
}

library("easypackages")
libraries(list.of.packages)
libraries(listOfBiocPackages)
