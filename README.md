# misaem

### Introduction

`misaem` is an implementation of methodology which performs statistical inference for logistic regression model with missing data. This method is based on likelihood, including:

1. Estimate the parameters of logistic regression by a stochastic approximation version of EM algorithm;
2. Estimation of parameters' variance based one Louis formula;
3. Model selection procedure based on BIC;
4. Prediction on a test set which may contain missing values.

### Installation of package 
Now you can install the package **misaem** from CRAN. 
```{r}
install.packages("misaem")
 ```
### Using the misaem package
Basicly,

1. `miss.saem` contains the procedure of estimation for parameters, as well as their variance, and observed likelihood.
2. `model_selection` aims at selecting a best model according to BIC.
3. `pred_saem` performs prediction on a test set which may contain missing values.

For more details, You can find the vignette, which illustrate the basic and further usage of misaem package:
```{r}
library(misaem)
vignette('misaem')
 ```

### Reference
Stochastic Approximation EM for Logistic regression with missing values (2018, Jiang W., Josse J., Lavielle M., Traumabase group)" [arxiv link](https://arxiv.org/abs/1805.04602).


