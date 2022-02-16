#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#
set -o nounset -o pipefail -o errexit
# set -o xtrace


# linear regression
today=`date +%Y-%m-%d`
random_number=$(shuf -i1-100000 -n1)

outputFile=""

probesets_or_gene_symbols=$1;
CSV_file=$2;
geoDatasetCode=$3;
featureName=$4;
firstCondition=$5;
secondCondition=$6;
outputFolder=$7;


echo -e "\n- - - - - - - - - - - - - - - - - - - - - - - - - ";
echo "Input arguments: "
echo "probesets_or_gene_symbols: $probesets_or_gene_symbols";
echo "CSV_file: $CSV_file";
echo "geoDatasetCode: $geoDatasetCode";
echo "featureName: $featureName";
echo "firstCondition: $firstCondition";
echo "secondCondition: $secondCondition";
echo "outputFolder: $outputFolder";
echo -e "- - - - - - - - - - - - - - - - - - - - - - - - - \n";

suggestedFolder=$HOME"/"$outputFolder;
mkdir -p $suggestedFolder
outputResultsFile=$suggestedFolder"/"$today"_rand"$random_number"_execution_results.txt"
outputLogFile=$suggestedFolder"/"$today"_rand"$random_number"_execution_log.txt"

echo "Callin' easyDifferentialGeneCoexpression now; the results will be saved in the "$suggestedFolder" folder"

Rscript easyDifferentialGeneCoexpressionInputParameters.r "$probesets_or_gene_symbols" "$CSV_file" "$geoDatasetCode" "$featureName" "$firstCondition" "$secondCondition" > $outputResultsFile 2> $outputLogFile;

echo "easyDifferentialGeneCoexpression -- The end"

#                                                                                                                                                                                                                                                                                                   
# # Rscript easyDifferentialGeneCoexpressionInputParameters.r "PROBESETS" "../data/dc_probeset_list03.csv" "GSE30201" "source_name_ch1" "Patient" "Normal"
# 
# Rscript easyDifferentialGeneCoexpressionInputParameters.r $studyFlag > $outputFile 2> $outputFile
