# easyDifferentialGeneCoexpression #

`easyDifferentialGeneCoexpression`: an easy method to compute  differential gene coexpression.

## Summary ##

A function that reads in the GEO code of a list of probesets or gene symbols, a gene expression dataset GEO accession code, the name of the dataset feature discriminating the two conditions for the differential coexpression, and the values of the two different conditions for the differential coexpression. If the input gene list is made of gene symbols, this package associates the probesets to these gene symbols, if found. Platforms available: GPL80, GPL8300, GPL80, GPL96, GPL570, GPL571, GPL20115, GPL1293, GPL6102, GPL6104, GPL6883, GPL6884, GPL13497, GPL14550, GPL17077, GPL6480. GEO: Gene Expression Omnibus. ID: identifier code. 
The GEO datasets are downloaded from the URL <https://ftp.ncbi.nlm.nih.gov/geo/series/>.

This function has been designed for beginners and users having limited experience with `R`.

## Installation ##

To run `easyDifferentialGeneCoexpression`, you need to have the following programs and packages installed in your computer:

* R (version > 4.0)
* R packages: `easypackages, geneExpressionFromGEO, dplyr, jetset`
* R Bioconductor packages `Biobase, org.Hs.eg.db, GEOquery, diffcoexp, annotate, sva`

You can install the `easyDifferentialGeneCoexpression` package and its dependencies from CRAN, and load it, with the following commands typed in the `R` terminal console:

    R
    install.packages("easyDifferentialGeneCoexpression", repos='http://cran.us.r-project.org')
    library("easyDifferentialGeneCoexpression")
    
If it is impossible to download the package and its dependencies from CRAN, you can download the `easyDifferentialGeneCoexpression` package from this GitHub repository, and then can execute the `installPackages.r` script that will install all the dependencies automatically:

    cd easyDifferentialGeneCoexpression
    R
    source("installPackages.r")
    
Afterwards,  you execute the `easyDifferentialGeneCoexpression.r` file from an R terminal:

    source("easyDifferentialGeneCoexpression.r")

## Execution instructions ##

To run `easyDifferentialGeneCoexpression`, a user needs a list of microarray probesets, the GEO accession code of the dataset, the name of the feature of the dataset discriminating the two conditions for the differential coexpress, the value of the first condition of that feature, and the value of the second condition of that feature.
The last parameter allows you to decide if you want the function to print messages during its operations or not.

## An example ##

We want to compute the differential coexpression of the probesets `200738_s_at, 217356_s_at, 206686_at, 226452_at, 223172_s_at, 223193_x_at, 224314_s_at, 230630_at`, and `202022_at` in the [GSE3268 gene expression dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE3268). Here are the commands we can use in a `R` shell environment:

    verboseFlag <- TRUE
    probesetList <- c("200738_s_at", "217356_s_at", "206686_at", "226452_at", "223172_s_at",
                                        "223193_x_at",  "224314_s_at", "230630_at", "202022_at")
    geoCode <- "GSE16237"
    easyDifferentialGeneCoexpression(probesetList,  geoCode, "outcome of the patient:ch1", "Died of disease", "Alive", verboseFlag)
    
## Contacts ##

The `easyDifferentialGeneCoexpression` package was developed by [Davide Chicco](https://www.DavideChicco.it). Questions should be
addressed to davidechicco(AT)davidechicco.it
