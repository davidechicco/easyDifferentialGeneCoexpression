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
listOfBiocPackages <- c("Biobase", "GEOquery", "diffcoexp", "annotate", "org.Hs.eg.db"  )

bioCpackagesNotInstalled <- which( !listOfBiocPackages %in% rownames(installed.packages()) )
if(length(bioCpackagesNotInstalled)) cat("package missing listOfBiocPackages[", bioCpackagesNotInstalled, "]: ", listOfBiocPackages[bioCpackagesNotInstalled], "\n", sep="")

# check there's still something left to install
if( length(bioCpackagesNotInstalled) ) {
    BiocManager::install(listOfBiocPackages[bioCpackagesNotInstalled])
}

library("easypackages")
libraries(list.of.packages)
libraries(listOfBiocPackages)

# "gene_assignment
platformsWithGene_assignmentField <- c("GPL11532", "GPL23126", "GPL6244")
                
# "Gene Symbol"
platformsWithGeneSpaceSymbolField <- c("GPL80", "GPL8300", "GPL80", "GPL96", "GPL570", "GPL571")
                
# "gene_symbol"
platformsWithGene_SymbolField <- c("GPL20115")
                
# "symbol"
platformsWithSymbolField <- c("GPL1293", "GPL6102", "GPL6104", "GPL6883", "GPL6884")
                 
# "GENE_SYMBOL
platformsWith_GENE_SYMBOL_Field <- c("GPL13497", "GPL14550", "GPL17077", "GPL6480")

SIGNIFICANCE_THRESHOLD <- 0.005


# fromProbesetToGeneSymbol()
fromProbesetToGeneSymbol <- function(thisProbeset, thisPlatform,  this_platform_ann_df, verbose=FALSE) {
    
    thisGeneSymbol <- NULL
    
    # thisGeneSymbol <- this_platform_ann_df[this_platform_ann_df$ID==thisProbeset, ]$"Gene Symbol"

    if(verbose == TRUE) cat("probeset ", thisProbeset, " for the microarray platform ", thisPlatform, "\n", sep="")
    
    if(thisGEOplatform %in% platformsWithGeneSpaceSymbolField) thisGeneSymbol <- this_platform_ann_df[this_platform_ann_df$"ID"==thisProbeset,]$"Gene Symbol"
    if(thisGEOplatform %in% platformsWithGene_SymbolField) thisGeneSymbol <- this_platform_ann_df[this_platform_ann_df$"ID"==thisProbeset,]$"gene_symbol"
    if(thisGEOplatform %in% platformsWithSymbolField) thisGeneSymbol <- this_platform_ann_df[this_platform_ann_df$"ID"==thisProbeset,]$"symbol"
    if(thisGEOplatform %in% platformsWith_GENE_SYMBOL_Field) thisGeneSymbol <- this_platform_ann_df[this_platform_ann_df$"ID"==thisProbeset,]$"GENE_SYMBOL"
    
     if(verbose == TRUE) cat("gene symbol found ", thisGeneSymbol, "\n", sep="")
     
     if(is.null(thisGeneSymbol) & verbose == TRUE) cat("no gene symbol found for", thisProbeset, "\n", sep="\t")
     
     return(thisGeneSymbol)
}


# main()

args <- commandArgs(trailingOnly=TRUE)

# csv_file_name <- "probeset_list01.csv"
# GSE_code <- "GSE16020"
# featureNameToDiscriminateConditions <- "characteristics_ch1"
# firstConditionName <- "control"
# secondConditionName <- "monocytopenia"

# Rscript coexpr_usage_example.r GENE_SYMBOLS "gene_symbols_list02.csv" "GSE16020"  "characteristics_ch1" "control" "monocytopenia"^C



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

probesets_flag <- grepl("probeset|Probeset|PROBESET", probesets_or_gene_symbols) %>% any()
gene_symbols_flag <- grepl("symbol|SYMBOL|GENE_SYMBOL|gene_symbol", probesets_or_gene_symbols) %>% any()

# gene expression download
gset <- GEOquery::getGEO(GSE_code,  GSEMatrix =TRUE, getGPL=FALSE)
	      
thisGEOplatform <- toString((gset)[[1]]@annotation)
if(length(gset) > 1) idx <- grep(thisGEOplatform, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

gset_expression <- gset%>% exprs()

gsetPhenoDataDF <- as(gset@phenoData, 'data.frame')

# random shuffle
gset_expression <- gset_expression[sample(nrow(gset_expression)),] 

thisGEOplatformJetSetCode <- NULL
if(thisGEOplatform=="GPL97" || (thisGEOplatform=="GPL96")) thisGEOplatformJetSetCode <- "hgu133a"
if(thisGEOplatform=="GPL570") thisGEOplatformJetSetCode <- "hgu133plus2"

# to implement: gene symbols file read and association of the 
list_of_probesets_to_select <- NULL

if(gene_symbols_flag == TRUE) {
    list_of_gene_symbols_to_select <- read.csv(csv_file_name, header=FALSE, sep=",", stringsAsFactors=FALSE)
    list_of_gene_symbols_to_select <-  as.vector(t(list_of_gene_symbols_to_select))
    
    cat("List of input gene symbols:\n")
    cat(list_of_gene_symbols_to_select, sep=", ")
    cat("\n")
    
    cat("Retrieving the probesets of the input gene symbols on the ", thisGEOplatform, " microarray platform\n", sep="")
    list_of_probesets_to_select_temp <- jmap(thisGEOplatformJetSetCode, symbol = list_of_gene_symbols_to_select)
    list_of_probesets_to_select_temp2 <- as.data.frame(list_of_probesets_to_select_temp)$list_of_probesets_to_select_temp
    list_of_probesets_to_select <- list_of_probesets_to_select_temp2[!is.na(list_of_probesets_to_select_temp2)]
    
    geneSymbolsWithoutProbesets_temp <- list_of_probesets_to_select_temp[is.na(list_of_probesets_to_select_temp)] %>% names()
    geneSymbolsWithoutProbesets <- toString(paste(geneSymbolsWithoutProbesets_temp, sep=" "))
    
    cat("The user inserted ", length(list_of_gene_symbols_to_select), " gene symbols\n", sep="")
    cat("The script will use ", length(list_of_probesets_to_select), " probesets (", geneSymbolsWithoutProbesets, " do not have a probeset on this platform)\n", sep="")
    
}

if(probesets_flag == TRUE) {
    list_of_probesets_to_select <- read.csv(csv_file_name, header=FALSE, sep=",", stringsAsFactors=FALSE)
    list_of_probesets_to_select <-  as.vector(t(list_of_probesets_to_select))
}

    probesets_flag <- !is.null(list_of_probesets_to_select)


cat("List of input probesets:\n")
cat(list_of_probesets_to_select, sep=", ")
cat("\n")



# healthy_controls_gene_expression <- gset_expression[, grepl("control", gset$"characteristics_ch1", fixed=TRUE)] 
# patients_gene_expression <- gset_expression[, grepl("monocytopenia", gset$"characteristics_ch1", fixed=TRUE)] 

first_condition_gene_expression <- gset_expression[, grepl(firstConditionName, gsetPhenoDataDF[, featureNameToDiscriminateConditions]
, fixed=TRUE)] 
second_condition_gene_expression <- gset_expression[, grepl(secondConditionName, gsetPhenoDataDF[, featureNameToDiscriminateConditions]
, fixed=TRUE)] 

numProbesets <- -1

probesets_flag <- TRUE

if(probesets_flag == TRUE) {

    numProbesets <- list_of_probesets_to_select %>% length()
    cat("We consider only ", numProbesets, " probesets indicated manually\n", sep="")
    coexpr_results <- coexpr(first_condition_gene_expression[list_of_probesets_to_select,], second_condition_gene_expression[list_of_probesets_to_select,], r.method = "pearson")

} 


significant_coexpressed_probeset_pairs <- coexpr_results[(order(coexpr_results$"p.diffcor") & coexpr_results$"p.diffcor" < SIGNIFICANCE_THRESHOLD),c("Gene.1", "Gene.2", "p.diffcor", "q.diffcor", "cor.diff")] %>% unique()
# %>% head()
 
rownames(significant_coexpressed_probeset_pairs) <- paste0(significant_coexpressed_probeset_pairs$"Gene.1", ",", significant_coexpressed_probeset_pairs$"Gene.2")

# cat("significantly differentially coexpressed gene pairs (threshold p-value < ", SIGNIFICANCE_THRESHOLD,") :\n", sep="")
# print(significant_coexpressed_probeset_pairs)

platform_ann <- annotate::readGEOAnn(GEOAccNum = thisGEOplatform)
platform_ann_df <- as.data.frame(platform_ann, stringsAsFactors=FALSE)

pb <- 1
significant_coexpressed_probeset_pairs$geneSymbolLeft <- ""
significant_coexpressed_probeset_pairs$geneSymbolRight <- ""
for(pb in 1:(significant_coexpressed_probeset_pairs %>% nrow()))
{
    significant_coexpressed_probeset_pairs[pb,]$geneSymbolLeft <- fromProbesetToGeneSymbol(significant_coexpressed_probeset_pairs[pb,]$"Gene.1", thisGEOplatform,  platform_ann_df)
    significant_coexpressed_probeset_pairs[pb,]$geneSymbolRight<- fromProbesetToGeneSymbol(significant_coexpressed_probeset_pairs[pb,]$"Gene.2", thisGEOplatform, platform_ann_df)
    
}

cat("\nTop coexpresseed pairs of genes based on cor.diff:\n")
print(significant_coexpressed_probeset_pairs[,c("geneSymbolLeft", "geneSymbolRight",  "p.diffcor")])
cat("\n")
