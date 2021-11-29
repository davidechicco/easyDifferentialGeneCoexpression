setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")
set.seed(11)

source("installPackages.r")


# main()

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
  stop("At least 4 argument must be supplied\n", call.=FALSE)
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

source("easyDifferentialGeneCoexpression.r")

gset00 <- GEOquery::getGEO(GSE_code,  GSEMatrix =TRUE, getGPL=FALSE)
platformCode <- toString((gset00)[[1]]@annotation)

list_of_probesets_to_select <-  probesetRetrieval(probesets_or_gene_symbols, csv_file_name, platformCode)

easyDifferentialGeneCoexpression(list_of_probesets_to_select, GSE_code, featureNameToDiscriminateConditions, firstConditionName, secondConditionName)
