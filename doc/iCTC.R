## ----first, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## ----setuplib-----------------------------------------------------------------
library(iCTC)

## ----data, message=FALSE,warning = FALSE,include=TRUE, cache=FALSE------------
cell_samples<-iCTC::raw_test_data$Clearcell_Polaris_sample_test
head(cell_samples)

## ----main, message=FALSE,warning = FALSE, include=TRUE, cache=FALSE-----------
results<-iCTC(cell_samples=cell_samples, cases = c(4,5,6))

## ----output, message=FALSE,warning = FALSE,include=TRUE, cache=FALSE----------
results$CTC_probabilistic_score

