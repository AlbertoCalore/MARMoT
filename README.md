
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MARMoT

<!-- badges: start -->
<!-- badges: end -->

The goal of MARMoT is to …

## Installation

You can install the development version of MARMoT like so:

``` r
# FILL THIS IN! HOW CAN PEOPLE INSTALL YOUR DEV PACKAGE?
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(MARMoT)
out = MARMoT(data = MARMoT_data, confounders = c("race", "age"), treatment = "hospital", n.cores = 1)
#> [1] "Absolute standardized bias - before balancing"
#>          Min 1st quartile       Median 3rd quartile          Max         Mean 
#>        0.010        2.875        6.859       12.027       50.700        8.667 
#>      Over 5%     Over 10% 
#>      214.000      118.000 
#> [1] "... Computing average rank ..."
#> [1] "... Balancing ..."
#> [1] "Absolute standardized bias - after balancing"
#>          Min 1st quartile       Median 3rd quartile          Max         Mean 
#>        0.132        0.138        0.412        0.570       43.302        1.055 
#>      Over 5%     Over 10% 
#>       23.000        4.000
```
