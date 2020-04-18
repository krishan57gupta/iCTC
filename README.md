
<!-- README.md is generated from README.Rmd. Please edit that file -->

## iCTC: identification of circulating tumor cells

The goal of iCTC is to detect whether peripheral blood cells have CTCs
(circulating tumor cell) or not.

## Installation

The developer version of the R package can be installed with the
following R commands:

``` r
library(devtools)
install_github("immunogenomics/harmony")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install('iCTC')
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
results<-iCTC(cell_samples=cell_samples, cases = c(4,5,6))
#> [1] "All desired genes in row names found"
#> [1] "All samples found with atlest 10 percent expressed genes\n                in column names of your data"
#> [1] 4 5 6
#> [1] "Harmony count 0"
#> [1] "PCA count 3"
#> [1] "Original count 0"
#> [1] "PCA correction running..."
#> [1] "Normalization Done"
#> [1] "Transformation Done"
#> [1] "PCA correction Done"
#> [1] "NB running..."
#> [1] "NB Done"
#> [1] "RF running..."
#> [1] "RF Done"
#> [1] "GBM running..."
#> [1] "GBM Done"
```

``` r
results$CTC_probabilistic_score
#>         Ctc_Naveen_1850-061-049-CS26_S19_15 Ctc_Naveen_1850-061-049-CS29_S22_15
#> PCA_NB                            1.0000000                           1.0000000
#> PCA_RF                            0.7100000                           0.7080000
#> PCA_GBM                           0.9997641                           0.9997536
#>         Ctc_Naveen_1851009049-10_S26_15 Ctc_Naveen_1851009049-11_S34_15
#> PCA_NB                         1.000000                       1.0000000
#> PCA_RF                         0.694000                       0.6820000
#> PCA_GBM                        0.999689                       0.9993584
#>         Ctc_Naveen_1851009049-14_S11_15 Ctc_Naveen_1851009049-15_S19_15
#> PCA_NB                        1.0000000                       1.0000000
#> PCA_RF                        0.6940000                       0.7080000
#> PCA_GBM                       0.9996159                       0.9995862
#>         Ctc_Naveen_1851009049-16_S27_15 Ctc_Naveen_1851009049-17_S35_15
#> PCA_NB                        1.0000000                       1.0000000
#> PCA_RF                        0.6740000                       0.6300000
#> PCA_GBM                       0.9993271                       0.9995179
#>         Ctc_Naveen_1851009049-4_S25_15 Ctc_Naveen_1851009049-5_S33_15
#> PCA_NB                       1.0000000                      1.0000000
#> PCA_RF                       0.6940000                      0.7100000
#> PCA_GBM                      0.9995179                      0.9995734
#>         Ctc_Naveen_1851009049-7_S2_15 Ctc_Naveen_1851009049-8_S10_15
#> PCA_NB                       1.000000                       1.000000
#> PCA_RF                       0.690000                       0.708000
#> PCA_GBM                      0.999336                       0.999767
#>         Ctc_Naveen_1851009049-9_S18_15 Ctc_Naveen_1851013039-26_S61_15
#> PCA_NB                       1.0000000                       1.0000000
#> PCA_RF                       0.6980000                       0.9220000
#> PCA_GBM                      0.9996666                       0.9725238
#>         Ctc_Naveen_1851013039-6_S89_15
#> PCA_NB                       1.0000000
#> PCA_RF                       0.5740000
#> PCA_GBM                      0.7853166
```
