# misaem
### Package R for "Stochastic Approximation EM for Logistic Regression with Missing Values (W. Jiang, J. Josse, M. Lavielle, Traumabase Group, 2018)"

`misaem` is a method to apply statistical inference for logistic regression model with missing data. This method is based on likelihood, including 
1. A stochastic approximation version of EM algorithm based on Metropolis-Hasting sampling, to estimate the parameters of logistic regression;
2. Estimation of parameters' variance based one Louis formula;
3. Model selection procedure based on BIC.

### Installation of package 
First you can install the package **misaem** from Github. The main function `miss.saem` contains the procedure of estimation for parameters, as well as their variance, and observed likelihood.
```{r}
library(devtools)
install_github("wjiang94/misaem")
 ```
### Using the misaem package
You can find the vignette, which illustrate the basic and further usage of misaem package:
```{r}
library(misaem)
vignette('misaem')
 ```

### Reference
Stochastic Approximation EM for Logistic regression with missing values (2018, Jiang W., Josse J., Lavielle M., Traumabase group)" [arxiv link](https://arxiv.org/abs/1805.04602).


