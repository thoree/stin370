
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Functions and data for STIN370

The purpose

## Installation

``` r
 # First install devtools if needed
if(!require(devtools)) install.packages("devtools")

# Install dvir from GitHub
devtools::install_github("thoree/stin370")
```

## Load libraries

``` r
library(stin370, quietly = T)
```

``` r
library(genetics, quietly = T) # May generate warnings
```

## Simulate data from unrelated individuals

``` r
nInd = 100
db = simulateUnrelated(nInd = nInd, nMark = 2, alleles = 1:2,
                       afreq = c(0.4, 0.6), seed = 1729)
```

## Estimate allele frequencies

Extract and show some genotypes from first marker

``` r
(res = estimateAlleleFrequencies(db, markers = 1))
#> $`Marker 1`
#>   Count Proportion
#> 2   129      0.645
#> 1    71      0.355
```

## HWE test

There are many options, e.g.,

``` r
g = getGenotypes(db)
g1 = genotype(g[,1])
(test1 = HWE.chisq(g1))
#> 
#>  Pearson's Chi-squared test with simulated p-value (based on 10000
#>  replicates)
#> 
#> data:  tab
#> X-squared = 1.2918, df = NA, p-value = 0.2834
```

Checking test-statistic

``` r
p1 = 0.355
p2 = 0.645
expected = nInd *c(p1^2, 2*p1*p2, p2^2)
observed = c(10, 51, 39)
df = data.frame(genotype = c("1/1", "1/2", "2/2"), observed, expected)
df
#>   genotype observed expected
#> 1      1/1       10  12.6025
#> 2      1/2       51  45.7950
#> 3      2/2       39  41.6025
```

``` r
sum((observed -expected)^2/expected)
#> [1] 1.29183
```
