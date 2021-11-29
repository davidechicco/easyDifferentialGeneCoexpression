setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")
set.seed(11)

probesetList <- c("200738_s_at", "17356_s_at", "206686_at", "226452_at", "223172_s_at",  "223193_x_at",  "224314_s_at", "230630_at", "202022_at")

easyDifferentialGeneCoexpression(probesetList,  "GSE16237", "outcome of the patient:ch1", "Died of disease", "Alive")

# Rscript main.r "GENE_SYMBOLS" "../data/gene_symbols_list02.csv" "GSE16020"  "characteristics_ch1" "control" "monocytopenia"

# Rscript main.r "PROBESETS" "../data/dc_probeset_list03.csv" "GSE16237" "outcome of the patient:ch1" "Died of disease" "Alive"
