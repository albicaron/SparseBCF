# Shrinkage Bayesian Causal Forests
This page gathers `R` code to implement "Shrinkage Bayesian Causal Forests for Heterogeneous Treatment Effects Estimation" as in [Caron et al. (2021)](https://arxiv.org/pdf/2102.06573.pdf). The `SparseBCF` webpage includes both the `R` package and the code to replicate simulated and real-world examples found in the paper, [Caron et al. (2021)](https://arxiv.org/pdf/2102.06573.pdf).

`SparseBCF` is a powerful Bayesian nonparametric regression model to estimate individualized/heterogeneous treatment effects, and design optimal treatment administration policies. It builds on top of the popular work on "Bayesian Causal Forests" by [Hahn, Murray, and Carvalho (2020)](https://projecteuclid.org/euclid.ba/1580461461), which has already been shown to perform extremely well on a variety of tasks. `SparseBCF` adds Dirichlet priors to the model in order to induce and adapt to sparsity, and most importantly to detect what are the main variables responsible for the heterogeneity behind treatment response. More details about the methodology can be found in the paper [Caron et al. (2021)](https://arxiv.org/pdf/2102.06573.pdf).

## Package Installation 
`SparseBCF` package requires compilation, so you must have "Rtools" installed - click [here](https://cran.r-project.org/bin/windows/Rtools/) for details on how to install it.

To install the package via Github:
```{r}
if (!require("devtools")) install.packages("devtools")

devtools::install_github("albicaron/SparseBCF")
```

## Examples 
The `Simulated examples` and `Real world example` folders contain code to replicate the examples displayed in the paper [Caron et al. (2021)](https://arxiv.org/pdf/2102.06573.pdf).
