
## iCTC: identification of circulating tumor cells

The goal of iCTC is to detect whether peripheral blood cells have CTCs
(circulating tumor cell) or not.

## Installation

The developer version of the R package can be installed with the
following R commands:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install('ROSeq')
```

or can be installed with the following R commands:

``` r
library(devtools)
install_github("immunogenomics/harmony")
install_github('krishan57gupta/iCTC')
```

## Vignette tutorial

This vignette uses a small dataset of cell samples, which is saved in
package itself, to demonstrate a standard pipeline. This vignette can be
used as a tutorial as well.

Libraries need to be loaded before running iCTC.

``` r
library(iCTC)
```

``` r
cell_samples=iCTC::raw_test_data$Clearcell_Polaris_sample_test
head(cell_samples)
#>          Ctc_Naveen_1850-061-049-CS26_S19_15
#> CDH2                                    0.00
#> SIGLEC14                                1.92
#> GNPDA1                                  0.00
#> KCNE3                                   3.81
#> CDH3                                    0.00
#> TANK                                    0.00
#>          Ctc_Naveen_1850-061-049-CS29_S22_15 Ctc_Naveen_1851009049-10_S26_15
#> CDH2                                    0.00                            0.00
#> SIGLEC14                                0.00                            1.00
#> GNPDA1                                  0.00                            0.00
#> KCNE3                                   2.78                            2.90
#> CDH3                                    0.00                            0.31
#> TANK                                    0.00                            0.00
#>          Ctc_Naveen_1851009049-11_S34_15 Ctc_Naveen_1851009049-14_S11_15
#> CDH2                                0.00                            0.00
#> SIGLEC14                            1.00                            1.70
#> GNPDA1                              0.00                            0.00
#> KCNE3                               6.09                            3.07
#> CDH3                                0.00                            0.00
#> TANK                                0.00                            0.00
#>          Ctc_Naveen_1851009049-15_S19_15 Ctc_Naveen_1851009049-16_S27_15
#> CDH2                                 0.0                            0.00
#> SIGLEC14                             0.0                            0.00
#> GNPDA1                               0.0                            0.00
#> KCNE3                                2.9                            5.68
#> CDH3                                 0.0                            0.00
#> TANK                                 0.0                            0.00
#>          Ctc_Naveen_1851009049-17_S35_15 Ctc_Naveen_1851009049-4_S25_15
#> CDH2                                0.00                           0.00
#> SIGLEC14                            1.00                           0.00
#> GNPDA1                              0.00                           0.00
#> KCNE3                              10.69                           2.52
#> CDH3                                0.00                           0.00
#> TANK                                0.00                           0.00
#>          Ctc_Naveen_1851009049-5_S33_15 Ctc_Naveen_1851009049-7_S2_15
#> CDH2                               0.00                           0.0
#> SIGLEC14                           3.07                           1.0
#> GNPDA1                             0.00                           0.0
#> KCNE3                              5.00                           5.8
#> CDH3                               0.00                           0.0
#> TANK                               0.00                           0.0
#>          Ctc_Naveen_1851009049-8_S10_15 Ctc_Naveen_1851009049-9_S18_15
#> CDH2                               0.00                           0.00
#> SIGLEC14                           0.00                           3.00
#> GNPDA1                             0.00                           0.00
#> KCNE3                              4.22                           4.01
#> CDH3                               0.00                           0.00
#> TANK                               0.00                           0.00
#>          Ctc_Naveen_1851013039-26_S61_15 Ctc_Naveen_1851013039-6_S89_15
#> CDH2                                0.00                           0.00
#> SIGLEC14                            0.00                           1.00
#> GNPDA1                              0.00                           0.00
#> KCNE3                               2.83                           0.97
#> CDH3                                0.00                           0.00
#> TANK                                0.00                           0.00
```

``` r
results<-iCTC(cell_samples=cell_samples, cases = seq_len(3))
#> [1] "All desired genes in row names found"
#> [1] "All samples found with atlest 10 percent expressed genes\n                in column names of your data"
#> [1] 1 2 3
#> [1] "Harmony count 3"
#> [1] "PCA count 0"
#> [1] "Original count 0"
#> [1] "Harmony correction running..."
#> [1] "Normalization Done"
#> [1] "Transformation Done"
#> [1] "Harmony correction Done"
#> [1] "NB running..."
#> [1] "NB Done"
#> [1] "RF running..."
#> [1] "RF Done"
#> [1] "GBM running..."
#> [1] "GBM Done"
```

``` r
results
#> $predicted_Labels
#>             Ctc_Naveen_1850-061-049-CS26_S19_15
#> Harmony_NB  "CTC"                              
#> Harmony_RF  "CTC"                              
#> Harmony_GBM "CTC"                              
#>             Ctc_Naveen_1850-061-049-CS29_S22_15 Ctc_Naveen_1851009049-10_S26_15
#> Harmony_NB  "CTC"                               "CTC"                          
#> Harmony_RF  "CTC"                               "CTC"                          
#> Harmony_GBM "CTC"                               "CTC"                          
#>             Ctc_Naveen_1851009049-11_S34_15 Ctc_Naveen_1851009049-14_S11_15
#> Harmony_NB  "CTC"                           "CTC"                          
#> Harmony_RF  "CTC"                           "CTC"                          
#> Harmony_GBM "CTC"                           "CTC"                          
#>             Ctc_Naveen_1851009049-15_S19_15 Ctc_Naveen_1851009049-16_S27_15
#> Harmony_NB  "CTC"                           "CTC"                          
#> Harmony_RF  "CTC"                           "CTC"                          
#> Harmony_GBM "CTC"                           "CTC"                          
#>             Ctc_Naveen_1851009049-17_S35_15 Ctc_Naveen_1851009049-4_S25_15
#> Harmony_NB  "CTC"                           "CTC"                         
#> Harmony_RF  "CTC"                           "CTC"                         
#> Harmony_GBM "CTC"                           "CTC"                         
#>             Ctc_Naveen_1851009049-5_S33_15 Ctc_Naveen_1851009049-7_S2_15
#> Harmony_NB  "CTC"                          "CTC"                        
#> Harmony_RF  "CTC"                          "CTC"                        
#> Harmony_GBM "CTC"                          "CTC"                        
#>             Ctc_Naveen_1851009049-8_S10_15 Ctc_Naveen_1851009049-9_S18_15
#> Harmony_NB  "CTC"                          "CTC"                         
#> Harmony_RF  "CTC"                          "CTC"                         
#> Harmony_GBM "CTC"                          "CTC"                         
#>             Ctc_Naveen_1851013039-26_S61_15 Ctc_Naveen_1851013039-6_S89_15
#> Harmony_NB  "CTC"                           "CTC"                         
#> Harmony_RF  "CTC"                           "Blood"                       
#> Harmony_GBM "CTC"                           "CTC"
```
