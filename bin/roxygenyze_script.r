
setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")
set.seed(11)


list.of.packages <- c("roxygen2") # other packages
new_packages_to_install <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new_packages_to_install)) install.packages(new_packages_to_install, repos="https://utstat.toronto.edu/cran/")

library("roxygen2")

setwd("../../easyDifferentialGeneCoexpression/")
cat(getwd(), "\n")

roxygen2::roxygenize()

cat("roxygenize() completed\n")
