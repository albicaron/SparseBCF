# Sparse Bayesian Causal Forests
This page gathers `R` code to implement "Sparse Bayesian Causal Forests for Heterogeneous Treatment Effects Estimation" as in CITE. This `SparseBCF` page includes both the `R` package and the simulated and real-world examples found in the paper CITE.

`SparseBCF` is a powerful Bayesian nonparametric regression model to estimate individualized treatment effects, and design optimal treatment administration policies. It builds on top of the popular work on "Bayesian Causal Forests" by [Hahn, Murray, and Carvalho 2020](https://projecteuclid.org/euclid.ba/1580461461), which has already been shown to perform very well on a variety of tasks, and adds Dirichlet priors to induce and adapt to sparsity - and most importantly to detect what are the main variables responsible for the heterogeneous treatment response. More details about the model can be found in the paper CITE.

## Package Installation 
`SparseBCF` package requires compilation, so you must have "Rtools" installed - click [here](https://cran.r-project.org/bin/windows/Rtools/) for details on how to install it.

To install the package via Github:
```{r}
if (!require("devtools")) install.packages("devtools")

devtools::install_github("albicaron/SparseBCF")
```

## Examples 
The `Simulated examples` and `Real world example` folders contain code to recreate the examples displayed in the paper CITE.
