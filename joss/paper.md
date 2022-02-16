---
title: '`easyDifferentialGeneCoexpression`, a handy bioinformatics tool to easily perform differential gene coexpression'
tags:
  - R
  - bioinformatics
  - computational biology
  - bioinformatics
  - gene expression
authors:
  - name: Davide Chicco
    orcid: 0000-0001-9655-7142
    affiliation: 1
  - name: Abbas Alameer
    orcid: 0000-0002-0699-163X
    affiliation: 2

affiliations:
 - name: University of Toronto
   index: 1
 - name: Kuwait University
   index: 2
date: 13 February 2022
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

Gene expression is a way to measure the activity of genes in a structured experiment: the more intense the gene expression is, the more active the gene is in that sample [@Chowdhury:2019].
Differential gene coexpression is a technique that indicates the pairs of genes that have different expression trends in samples having two different conditions, which then result in different correlation coefficients [@zheng2014gene].

If two genes have very different gene expression trends in data samples of patients with breast cancer and in data samples of healthy controls, it is possible that both those genes (or at least one of them) might have a significant role in breast cancer development and/or prognosis [@gov2017differential; @choi2005differential].

In this short paper, we present `easyDifferentialGeneCoexpression`, an R package which handily computes the differential gene expression between genes in a specific dataset, and returns the list of significant differential pairs of genes, if found.
Our software is available as a [R library](https://metacpan.org/pod/App::easyDifferentialGeneCoexpression), as a [Perl application](https://metacpan.org/pod/App::easyDifferentialGeneCoexpressionWrapper) that can be used in any standard terminal shell, and as a [repository on GitHub](https://github.com/davidechicco/easyDifferentialGeneCoexpression).

# Statement of need

Several R packages for differential gene expression already exist: `diffcoexp` [@WWW-diffcoexp; @yang2013dcgl], `decode` [@lui2015decode] and `dcanr` [@bhuva2019differential] on Bioconductor [@Huber2015]. 
However, they all have limitations, that can make differential gene coexpression analyses difficult, especially for beginners.
First, they provide results measured with multiple coefficients, which can be an asset for experienced researchers, but can also be confusing for beginners and unexperienced users.

The `diffcoexp()` function of the `diffcoexp` [@WWW-diffcoexp] library, for example, returns the pairs of differentially coexpressed genes ranked by difference between correlation coefficients under the second condition and the first condition (`cor.diff`), *p*-value under null hypothesis that difference between two correlation coefficients under two conditions equals to zero using Fisher’s r-to-Z transformation (`p.diffcor`), and adjusted *p*-value under null hypothesis that difference between two correlation coefficients under two conditions equals to zero using Fisher’s r-to-Z transformation (`q.diffcor`).
These three coefficients have different meanings and can generate three different rankings.
 Instead, our `easyDifferentialGeneCoexpression` package, that we built right on `diffcoexp`, generates a final ranking of pairs of significantly expressed genes only through the *p*-value difference ranking, which we believe it is the most informative coefficent and ranking.
 
 Additionally, our `easyDifferentialGeneCoexpression` package returns a list of significantly coexpressed gene pairs only if their *p*-values are stricly lower than the 0.005 significance threshold, as suggested by @benjamin2018redefine.
 To avoid *p*-hacking [@head2015extent], the users cannot choose their  prefered significance threshold. 
 By using this $5 \times 10^{-3}$ threshold, in fact, users can rest assured that any pair of coexpressed genes is significant enough to be reliable, avoiding insignificant results that could lead to unimportant discoveries [@ioannidis2005most].

 
# Example

To install `easyDifferentialGeneCoexpression` from CRAN, in an R environment:

```
    install.packages("easyDifferentialGeneCoexpression")
```

To install `easyDifferentialGeneCoexpression` from GitHub:

```
    git clone https://github.com/davidechicco/easyDifferentialGeneCoexpression.git
````

To install `easyDifferentialGeneCoexpression` from CPAN, on a Linux operating system:

```
    cpanm App::easyDifferentialGeneCoexpressionWrapper
```

Please notice that in a Linux Ubuntu system the user might have to run the last command in the `sudo` mode.

To use `easyDifferentialGeneCoexpression` in an R environment:

```
    ## Load the library
    library("easyDifferentialGeneCoexpression")
    
    ## List of probesets of the genes for which
    ## to compute the differential gene expression
    probesetList <- c("200738_s_at", "217356_s_at",
    "206686_at", "226452_at", "223172_s_at", "223193_x_at", 
    "224314_s_at", "230630_at", "202022_at")
    
    ## Function parameters
    verboseFlag <- TRUE
    datasetGEOcode <- "GSE16237"
    conditionFeatureName <- "outcome of the patient:ch1"
    firstConditionName <-  "Died of disease"
    secondConditionName <- "Alive"
    batchCorrectionFlag <-  TRUE
    
    ## Function call
    easyDifferentialGeneCoexpression(probesetList, datasetGEOcode, conditionFeatureName, 
    firstConditionName, secondConditionName, batchCorrectionFlag, verboseFlag)
```

The output of the call is the following result:

```
    Significant top coexpresseed pairs of genes based
    on p-value difference (p.diffcor < 0.005):

                            geneSymbolLeft geneSymbolRight p.diffcor
    223172_s_at,206686_at            MTFP1 PDK1         3.111242e-06
    223193_x_at,206686_at          FAM162A PDK1         1.022005e-05
    223193_x_at,223172_s_at        FAM162A MTFP1        1.132584e-05
    217356_s_at,223172_s_at           PGK1 MTFP1        3.917640e-05
    226452_at,206686_at               PDK1 PDK1         1.956924e-04
    217356_s_at,223193_x_at           PGK1 FAM162A      7.650626e-04
    217356_s_at,226452_at             PGK1 PDK1         1.132578e-03
    226452_at,223172_s_at             PDK1 MTFP1        2.243446e-03
    200738_s_at,226452_at             PGK1 PDK1         2.506663e-03
    217356_s_at,206686_at             PGK1 PDK1         4.742024e-03


```

In this example, we computed the differential gene coexpression of the genes related to the probesets saved in the `probesetList` variable in the [GSE16237 gene expression dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16237) [@ohtaki2010robust] available on Gene Expression Omnibus (GEO).
Other prognostic  gene expression datasets for this scope can be found through the recently released Perl package `geoCancerPrognosticDatasetsRetriever` [@alameer2021geocancerprognosticdatasetsretriever].

`easyDifferentialGeneCoexpression` accepts probesets or gene symbols as input; in the latter case, it associates the input gene symbols to the corresponding microarray platform probesets through the `annotate` [@WWW-PACKAGE-ANNOTATE] and the `geneExpressionFromGEO` packages [@chicco2022geneexpressionfromgeo].

The GSE16237 dataset contains prognostic gene expression samples of 51 patients diagnosed with neuroblastoma. 
In this cohort, 39 patients died of this childhood cancer and 12 patients survived. This condition is encoded in the `"outcome of the patient:ch1"` variable of the dataset: the `"Died of disease"` label indicates the deceased patients and the `"Alive"` label indicates the survived individuals, of course.
In the reported R code example, we specified all these pieces of information in the `datasetGEOcode`, `conditionFeatureName`, `firstConditionName`, and `secondConditionName` variables.

Our package main function `easyDifferentialGeneCoexpression()` first performs a batch correction [@chen2011removing],if the batch correction paramenter is set to true, and then computes the differential gene coexpression. This function eventually generates a list of significantly coexpressed pairs of genes, whose *p*-value is lower than 0.005, as suggested by @benjamin2018redefine.
To avoid *p*-hacking [@head2015extent], this threshold cannot be changed by the user.

In the results, the MTFP1-PDK1, FAM162A-PDK1, and MTFP1-FAM162A gene pairs result being the most significantly coexpressed gene pairs, suggesting an active role of these three genes (FAM162A, MTFP1, and PDK1) in neuroblastoma. Researchers can use this information to carry on new experiments and scientific analyses investigating the role of these three genes in neuroblastoma.

# Acknowledgements

The authors thank the CPAN, and CRAN community members for their help.

# References
