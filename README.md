# jeffreysmeta

<!-- badges: start -->
[![R-CMD-check](https://github.com/mathurlabstanford/jeffreysmeta/workflows/R-CMD-check/badge.svg)](https://github.com/mathurlabstanford/jeffreysmeta/actions)
<!-- badges: end -->
  
`jeffreysmeta` is an R package that provides a bias correction for the joint effects of p-hacking (i.e., manipulation of results within studies to obtain significant, positive estimates) and traditional publication bias (i.e., the selective publication of studies with significant, positive results) in meta-analyses (per [Mathur, 2022](https://osf.io/ezjsx/)).

## Installation

You can install jeffreysmeta from CRAN with:
```
install.packages("jeffreysmeta")
```

You can install the development version of jeffreysmeta from [GitHub](https://github.com/) with:
``` r
# install.packages("devtools")
devtools::install_github("mathurlabstanford/jeffreysmeta")
```

You may also need to install [Stan](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).

## Example

Fit a bias-corrected meta-analysis of an example dataset from the package.

``` r
library(jeffreysmeta)
jeffreysmeta_meta(money_priming_meta$yi, money_priming_meta$vi)
```
