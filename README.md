
# prepost

<!-- badges: start -->

[![R-CMD-check](https://github.com/mattblackwell/prepost/workflows/R-CMD-check/badge.svg)](https://github.com/mattblackwell/prepost/actions)
<!-- badges: end -->

## Overview

This package contains code and sample data to implement the
non-parametric bounds and Bayesian methods for assessing priming and
post-treatment bias in experimental studies under various assumptions.

To get started, please see the article that developed these methods:

- Blackwell, Matthew et al (2025). Priming Bias Versus Post-Treatment
  Bias in Experimental Designs. *Political Analysis* Published online
  2025:1-17. <doi:10.1017/pan.2025.3>. ([journal
  version](http://doi.org/10.1017/pan.2025.3), [arXiv
  preprint](https://arxiv.org/pdf/2306.01211))

## Installation

``` r
## Install developer version
## install.packages("devtools")
devtools::install_github("mattblackwell/prepost", build_vignettes = TRUE)
```

## Usage

Both the nonparametric and Bayesian estimators all have prefixes that
indicate what type of experimental design being used.

- `pre_` functions can analyze data from a **pre-test design** where the
  moderator is measured pre-treatment.
- `post_` functions can analyze data from a **post-test design** where
  the moderator is measured post-treatment.
- `prepost_` functions can analyze data from a **random placement
  design**, in which the moderator is randomly assigned to be measured
  before or after treatment.

Most functions can be specified with a formula to identify the outcome
and treatment and another one-sided formula for the moderator:

``` r
library(prepost)
data(delponte)
out <- pre_bounds(
  formula = angry_bin ~ t_commonality,
   data = delponte,
  moderator = ~ itaid_bin
)
out
```

    ## $lower
    ##            
    ## -0.5923203 
    ## 
    ## $upper
    ##           
    ## 0.3221525 
    ## 
    ## $ci_lower
    ## [1] -0.6875343
    ## 
    ## $ci_upper
    ## [1] 0.4035053
    ## 
    ## $pre_est
    ## [1] -0.2701678
