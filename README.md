
<!-- README.md is generated from README.Rmd. Please edit that file -->

# METER

<!-- badges: start -->
<!-- badges: end -->

METER is an R package designed for analyzing tumor content (TC) in
cell-free DNA (cfDNA) samples using low-pass whole-genome bisulfite
sequencing (lpWGBS) data. It leverages differentially methylated sites
(DMS) and regions (DMR) to provide a comprehensive computational
framework. The package includes three modules: \* **METER-quant** to
measure ctDNA, based on tumor-specific DMS \* **METER-detect** to detect
ctDNA, based on tumor-specific DMR \* **METER-subtype** to infer tumor
subtype from ctDNA, based on tumor subtype-specific DMR

## Installation

You can install the development version of METER from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("caos-lab-unifi/METER")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(METER)
## basic example code
```

<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->
<!-- ```{r cars} -->
<!-- summary(cars) -->
<!-- ``` -->
<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. -->
<!-- You can also embed plots, for example: -->
<!-- ```{r pressure, echo = FALSE} -->
<!-- plot(pressure) -->
<!-- ``` -->
<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
