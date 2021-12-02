setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")
set.seed(11)



args <- commandArgs(trailingOnly=TRUE)

# csv_file_name <- "probeset_list01.csv"
# GSE_code <- "GSE16020"
# featureNameToDiscriminateConditions <- "characteristics_ch1"
# firstConditionName <- "control"
# secondConditionName <- "monocytopenia"

# Rscript main.r GENE_SYMBOLS "gene_symbols_list02.csv" "GSE16020"  "characteristics_ch1" "control" "monocytopenia"

probesets_or_gene_symbols <- ""
csv_file_name <- ""
GSE_code <- ""
featureNameToDiscriminateConditions <- ""
firstConditionName <- ""
secondConditionName <- ""


# test if there is at least one argument: if not, return an error
if (length(args) != 6) {
  stop("To run this package, 6 argument must be supplied\n", call.=FALSE)
} else {

  probesets_or_gene_symbols <-  toString(args[1])
  csv_file_name <- toString(args[2])
  GSE_code <- toString(args[3])
  featureNameToDiscriminateConditions <- toString(args[4])
  firstConditionName <- toString(args[5])
  secondConditionName <- toString(args[6])
  
}

cat("- - - - - - - - - - - - - - - - - - - - - - - - - \n")
cat("Input arguments values:\n")
cat("probesets_or_gene_symbols: ", probesets_or_gene_symbols, "\n", sep="")
cat("csv_file_name: ", csv_file_name, "\n", sep="")
cat("GSE_code: ", GSE_code, "\n", sep="")
cat("featureNameToDiscriminateConditions: ", featureNameToDiscriminateConditions, "\n", sep="")
cat("firstConditionName: ", firstConditionName, "\n", sep="")
cat("secondConditionName: ", secondConditionName, "\n", sep="")
cat("- - - - - - - - - - - - - - - - - - - - - - - - - \n")

probesets_or_gene_symbols_original <- probesets_or_gene_symbols
probesets_or_gene_symbols <- tolower(probesets_or_gene_symbols)

accepted_keywords <- c("probeset", "probesets", "probe_set",  "probe_sets", "probe sets",  "probe set", "symbol", "symbols", "gene symbol", "gene_symbol", "gene symbols", "gene_symbols")

if(!(probesets_or_gene_symbols %in% accepted_keywords)){

    cat("Error: the input probeset / gene symbols argument is: ", probesets_or_gene_symbols, "\n", sep="")
    cat("The input probeset / gene symbols argument should be one of the following terms (uppercase or lowercase): ", sep="")
    cat(accepted_keywords, sep=", ")
    cat("\nThe program will terminate here\n")
    quit(save="no")

}

source("installPackages.r")
source("easyDifferentialGeneCoexpression.r")

# list.of.packages <- c("easyDifferentialGeneCoexpression") # other packages
# new_packages_to_install <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new_packages_to_install)) install.packages(new_packages_to_install, repos="https://utstat.toronto.edu/cran/")
# 
# library("easyDifferentialGeneCoexpression")

gset00 <- geoDataDownload(GSE_code)
platformCode <- toString((gset00)[[1]]@annotation)
verboseFlag <- TRUE

list_of_probesets_to_select <-  probesetRetrieval(probesets_or_gene_symbols, csv_file_name, platformCode, verboseFlag)

easyDifferentialGeneCoexpressionResults <-  easyDifferentialGeneCoexpression(list_of_probesets_to_select, GSE_code, featureNameToDiscriminateConditions, firstConditionName, secondConditionName, verboseFlag)
