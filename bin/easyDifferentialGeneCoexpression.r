
#' Function that associates a gene symbol to a probeset for some Affymetrix platforms 
#'
#' @param thisProbeset probeset in input
#' @param thisPlatform GEO platform accession code
#' @param this_platform_ann_df annotation dataframe of the platform
#' @param verbose annotation verbose flag
#' @export
#' @import GEOquery
#' @return a gene symbol as string 
fromProbesetToGeneSymbol <- function(thisProbeset, thisPlatform,  this_platform_ann_df, verbose=FALSE) {
    
    thisGeneSymbol <- NULL
    
    platformsWithGeneSpaceSymbolField <- c("GPL80", "GPL8300", "GPL80", "GPL96", "GPL570", "GPL571") # "Gene Symbol"
    platformsWithGene_SymbolField <- c("GPL20115") # "gene_symbol"
    platformsWithSymbolField <- c("GPL1293", "GPL6102", "GPL6104", "GPL6883", "GPL6884") # "symbol"
    platformsWith_GENE_SYMBOL_Field <- c("GPL13497", "GPL14550", "GPL17077", "GPL6480") # "GENE_SYMBOL
    
    if(!(thisPlatform %in% c(platformsWithGeneSpaceSymbolField, platformsWithGene_SymbolField, platformsWithSymbolField, platformsWith_GENE_SYMBOL_Field))) {
    
        cat("The input platform ", thisPlatform, " is not among the ones available, the probeset gene symbol mapping is impossible.\n", sep="")
        return(thisGeneSymbol)
    }
    
    # thisGeneSymbol <- this_platform_ann_df[this_platform_ann_df$ID==thisProbeset, ]$"Gene Symbol"

    if(verbose == TRUE) cat("probeset ", thisProbeset, " for the microarray platform ", thisPlatform, "\n", sep="")
    
    if(thisPlatform %in% platformsWithGeneSpaceSymbolField) thisGeneSymbol <- this_platform_ann_df[this_platform_ann_df$"ID"==thisProbeset,]$"Gene Symbol"
    else if(thisPlatform %in% platformsWithGene_SymbolField) thisGeneSymbol <- this_platform_ann_df[this_platform_ann_df$"ID"==thisProbeset,]$"gene_symbol"
    else if(thisPlatform %in% platformsWithSymbolField) thisGeneSymbol <- this_platform_ann_df[this_platform_ann_df$"ID"==thisProbeset,]$"symbol"
    else if(thisPlatform %in% platformsWith_GENE_SYMBOL_Field) thisGeneSymbol <- this_platform_ann_df[this_platform_ann_df$"ID"==thisProbeset,]$"GENE_SYMBOL"
    
     if(verbose == TRUE) cat("gene symbol found ", thisGeneSymbol, "\n", sep="")
     
     if(is.null(thisGeneSymbol) & verbose == TRUE) cat("no gene symbol found for", thisProbeset, "\n", sep="\t")
     
     return(thisGeneSymbol)
}

#' Function that reads a CSV file of probesets or gene symbols and, in the latter case, it retrieves the original probesets
#'
#' @param probesets_or_gene_symbols flag saying if we're reading probesets or gene symbols
#' @param csv_file_name complete name of CSV file containing the probesets or the gene symbols
#' @param GEOcode
#' @export
#' @import GEOquery
#' @return a vector of probesets
probesetRetrieval <- function(probesets_or_gene_symbols, csv_file_name, platformCode) {

        probesets_flag <- grepl("probeset|Probeset|PROBESET", probesets_or_gene_symbols) %>% any()
        gene_symbols_flag <- grepl("symbol|SYMBOL|GENE_SYMBOL|gene_symbol", probesets_or_gene_symbols) %>% any()

        
        thisGEOplatformJetSetCode <- NULL
        if(platformCode=="GPL97" || (platformCode=="GPL96")) thisGEOplatformJetSetCode <- "hgu133a"
        else if(platformCode=="GPL570") thisGEOplatformJetSetCode <- "hgu133plus2"
        else { 
                cat("The platform of this dataset is not among the ones listed by Jetset. The program will stop here.") 
                quit(save="no")

            }

        # to implement: gene symbols file read and association of the 
        list_of_probesets_to_select <- NULL

        if(gene_symbols_flag == TRUE) {
            list_of_gene_symbols_to_select <- read.csv(csv_file_name, header=FALSE, sep=",", stringsAsFactors=FALSE)
            list_of_gene_symbols_to_select <-  as.vector(t(list_of_gene_symbols_to_select))
            
            cat("List of input gene symbols:\n")
            cat(list_of_gene_symbols_to_select, sep=", ")
            cat("\n")
            
            cat("Retrieving the probesets of the input gene symbols on the ", platformCode, " microarray platform\n", sep="")
            list_of_probesets_to_select_temp <- jmap(thisGEOplatformJetSetCode, symbol = list_of_gene_symbols_to_select)
            list_of_probesets_to_select_temp2 <- as.data.frame(list_of_probesets_to_select_temp)$list_of_probesets_to_select_temp
            list_of_probesets_to_select <- list_of_probesets_to_select_temp2[!is.na(list_of_probesets_to_select_temp2)]
            
            geneSymbolsWithoutProbesets_temp <- list_of_probesets_to_select_temp[is.na(list_of_probesets_to_select_temp)] %>% names()
            geneSymbolsWithoutProbesets <- toString(paste(geneSymbolsWithoutProbesets_temp, sep=" "))
            
            cat("The user inserted ", length(list_of_gene_symbols_to_select), " gene symbols\n", sep="")
            cat("The script will use ", length(list_of_probesets_to_select), " probesets (", geneSymbolsWithoutProbesets, " do not have a probeset on this platform)\n", sep="")
            
        } else if(probesets_flag == TRUE) {
            list_of_probesets_to_select <- read.csv(csv_file_name, header=FALSE, sep=",", stringsAsFactors=FALSE)
            list_of_probesets_to_select <-  as.vector(t(list_of_probesets_to_select))
        }

        probesets_flag <- !is.null(list_of_probesets_to_select)


        cat("List of input probesets:\n")
        cat(list_of_probesets_to_select, sep=", ")
        cat("\n")
        
        return(list_of_probesets_to_select)
        
        
}

#' Function that computes the differential coexpression of a list of probesets in a specific dataset and returns the most significant pairs
#'
#' @param list_of_probesets_to_select list of probesets for which the differential coexpression should be computed
#' @param GSE_code GEO accession code of the dataset to analyze
#' @param featureNameToDiscriminateConditions name of the feature of the dataset that contains the two conditions to investigate
#' @param firstConditionName name of the first condition in the feature to discriminate (for example, "healthy")
#' @param secondConditionName name of the second condition in the feature to discriminate (for example, "cancer")
#' @export
#' @import GEOquery annotate
#' @return a dataframe containing the significantly differentially co-expressed pairs of genes
#' @examples
#' probesetList <- c("242116_x_at", "225917_at", "227481_at", "210772_at", "218603_at", "225793_at", "212566_at", 
#' "43544_at", "218809_at", "202875_s_at", "233341_s_at", "228900_at")
#' signDiffCoexpressGenePairs <- easyDifferentialGeneCoexpression(probesetList,  "GSE16020"  "characteristics_ch1" "control" "monocytopenia")
easyDifferentialGeneCoexpression <- function(list_of_probesets_to_select, GSE_code, featureNameToDiscriminateConditions, firstConditionName, secondConditionName) 
{

        SIGNIFICANCE_THRESHOLD <- 0.005

        # gene expression download
        gset <- GEOquery::getGEO(GSE_code,  GSEMatrix =TRUE, getGPL=FALSE)
        thisGEOplatform <- toString((gset)[[1]]@annotation)
        
        if(length(gset) > 1) idx <- grep(thisGEOplatform, attr(gset, "names")) else idx <- 1
        gset <- gset[[idx]]

        gset_expression <- gset%>% exprs()
        gsetPhenoDataDF <- as(gset@phenoData, 'data.frame')

        # random shuffle
        gset_expression <- gset_expression[sample(nrow(gset_expression)),] 

        # healthy_controls_gene_expression <- gset_expression[, grepl("control", gset$"characteristics_ch1", fixed=TRUE)] 
        # patients_gene_expression <- gset_expression[, grepl("monocytopenia", gset$"characteristics_ch1", fixed=TRUE)] 

        first_condition_gene_expression <- gset_expression[, grepl(firstConditionName, gsetPhenoDataDF[, featureNameToDiscriminateConditions]
        , fixed=TRUE)] 
        second_condition_gene_expression <- gset_expression[, grepl(secondConditionName, gsetPhenoDataDF[, featureNameToDiscriminateConditions]
        , fixed=TRUE)] 
        
        cat("first_condition_gene_expression %>% dim()\n")
        print(first_condition_gene_expression %>% dim())
        cat("second_condition_gene_expression %>% dim()\n")
        print(second_condition_gene_expression %>% dim())

        numProbesets <- -1
        
        sharedProbesets <- intersect(rownames(gset_expression), list_of_probesets_to_select)
        unsharedProbesets <-  setdiff(list_of_probesets_to_select, rownames(gset_expression))

        numProbesets <- sharedProbesets %>% length()
        if(unsharedProbesets %>% length() >= 1) cat("Only ")
        cat(numProbesets, " of the ", list_of_probesets_to_select %>% length() ," input probesets are present in this dataset\n", sep="")
        if(unsharedProbesets %>% length() >= 1)  { 
                cat("The absent probesets are ", unsharedProbesets %>% length(),": ", sep="") 
                print(unsharedProbesets)
            }
        
        coexpr_results <- coexpr(first_condition_gene_expression[sharedProbesets,], second_condition_gene_expression[sharedProbesets,], r.method = "pearson")

        cat("Coexpression significance threshold: ", SIGNIFICANCE_THRESHOLD, "\n", sep="")
        
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
        for(pb in 1:(significant_coexpressed_probeset_pairs %>% nrow()))    {
            significant_coexpressed_probeset_pairs[pb,]$geneSymbolLeft <- fromProbesetToGeneSymbol(significant_coexpressed_probeset_pairs[pb,]$"Gene.1", thisGEOplatform,  platform_ann_df, TRUE)
            significant_coexpressed_probeset_pairs[pb,]$geneSymbolRight<- fromProbesetToGeneSymbol(significant_coexpressed_probeset_pairs[pb,]$"Gene.2", thisGEOplatform, platform_ann_df, TRUE)
            
        }
        
        colnames(significant_coexpressed_probeset_pairs)[1] <- c("probesetLeft")
        colnames(significant_coexpressed_probeset_pairs)[2] <- c("probesetRight")
        
        cat("\nTop coexpresseed pairs of genes based on cor.diff:\n")
        print(significant_coexpressed_probeset_pairs[,c("geneSymbolLeft", "geneSymbolRight",  "p.diffcor")])
        cat("\n\n")
        
        return(significant_coexpressed_probeset_pairs)

}
