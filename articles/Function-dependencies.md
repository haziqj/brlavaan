# Function dependencies

A useful graph to visualize the dependencies of the
[brlavaan](https://github.com/haziqj/sem-bias) package.

``` r
library(brlavaan)
#> Loading required package: lavaan
#> This is lavaan 0.6-21
#> lavaan is FREE software! Please report any bugs.
library(pkgdepR)
v <- deps("brlavaan")
print(v)
#> 
#> pkgdepR object
#> ------------------------------
#> Packages:    brlavaan
#> Total nodes: 71
#> Total links: 84
#>   -Between packages: 0
#>   -Within packages:  84
#>     --Between functions: 83
#>     --Self-referential:  1
plot(v)
```
