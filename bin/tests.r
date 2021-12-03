setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")
set.seed(11)

source("installPackages.r")
source("easyDifferentialGeneCoexpression.r")

probesetList <- c("200738_s_at", "217356_s_at", "206686_at", "226452_at", "223172_s_at",  "223193_x_at",  "224314_s_at", "230630_at", "202022_at")

# easyDifferentialGeneCoexpression(probesetList,  "GSE16237", "outcome of the patient:ch1", "Died of disease", "Alive")

# Rscript easyDifferentialGeneCoexpressionInputParameters.r  "GENE_SYMBOLS" "../data/gene_symbols_list02.csv" "GSE16020"  "characteristics_ch1" "control" "monocytopenia"

# Rscript easyDifferentialGeneCoexpressionInputParameters.r "PROBESETS" "../data/probeset_list01.csv "GSE16237" "outcome of the patient:ch1" "Died of disease" "Alive"

# Rscript easyDifferentialGeneCoexpressionInputParameters.r "PROBESETS" "../data/dc_probeset_list03.csv" "GSE30201" "source_name_ch1" "Patient" "Normal"

# easyDifferentialGeneCoexpression(probesetList,  "GSE30201", "source_name_ch1", "Patient", "Normal")


easyDifferentialGeneCoexpression(probesetList,  "GSE3268", "description", "Normal", "Tumor")

