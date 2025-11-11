
<!-- README.md is generated from README.Rmd. Please edit that file -->

# myregression

*A lightweight teaching package illustrating least squares estimation
via optimization and bootstrap inference.*

------------------------------------------------------------------------

## Overview

`myregression` provides a simple example of fitting linear models using
**numerical optimization (`optim`)** instead of the analytical OLS
formula, and demonstrates how to obtain **bootstrap standard errors**
and confidence intervals.

The package defines:

- `ls_loss()` — least-squares loss function  
- `my_lm()` — estimator using `optim()` and bootstrap covariance  
- S3 methods: `print()`, `coef()`, `vcov()`, `fitted()`, `residuals()`,
  `predict()`, `summary()`, `plot_marginals()`

This is ideal for illustrating numerical optimization, resampling, and
custom S3 interfaces in R.

------------------------------------------------------------------------

## Installation

``` r
# Install development tools if needed
install.packages("devtools")

# Install from GitHub
devtools::install_github("R-4-Data-Science/myregression")
```

------------------------------------------------------------------------

## Quick Example

``` r
library(myregression)

set.seed(1)
n <- 100; p <- 3
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
beta_true <- c(2, -1, 0.5)
y <- as.vector(cbind(1, X) %*% c(1, beta_true) + rnorm(n))

fit <- my_lm(X, y, B = 50)
print(fit)
#> Call:
#> my_lm(X = X, y = y, B = 50)
#> 
#> Coefficients:
#> (Intercept)          x1          x2          x3 
#>   1.0527311   1.9420411  -1.0549419   0.6046231 
#> 
#> Converged: TRUE  |  Bootstrap reps (kept): 50
summary(fit)
#> Call:
#> my_lm(X = X, y = y, B = 50)
#> 
#> Bootstrap-based coefficient table:
#>              Estimate Std. Error Statistic  Pr(>|t|)    
#> (Intercept)  1.052731   0.101006   10.4224 < 2.2e-16 ***
#> x1           1.942041   0.146496   13.2566 < 2.2e-16 ***
#> x2          -1.054942   0.116381   -9.0646 1.521e-14 ***
#> x3           0.604623   0.096831    6.2441 1.158e-08 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Confidence intervals (95%):
#>             CI Lower CI Upper
#> (Intercept)  0.85224  1.25323
#> x1           1.65125  2.23283
#> x2          -1.28596 -0.82393
#> x3           0.41241  0.79683
#> 
#> Residual df: 96  | Bootstrap reps used: 50
```

------------------------------------------------------------------------

## Plotting Marginals

``` r
plot_marginals(fit)
```

Each panel shows the scatter of `y` versus a predictor `x_j` along with
the *marginal fitted line* that holds other variables at their sample
means.

------------------------------------------------------------------------

## License

MIT © 2025 Roberto Molinari

------------------------------------------------------------------------

## Acknowledgments

Developed for instructional use in statistical computing and regression
analysis courses at **Auburn University**.
