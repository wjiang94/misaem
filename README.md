# misaem package

### Introduction

`misaem` is a package to perform linear regression and logistic regression with missing data, under MCAR (Missing completely at random) and MAR (Missing at random) mechanisms. The covariates are assumed to be continuous variables. The methodology implemented is based on maximization of the observed likelihood using EM-types of algorithms. The package includes:

1. Parameters estimation.
2. Estimation of standard deviation for estimated parameters.
3. Model selection procedure based on BIC. 

### Installation of package 
Now you can install the package **misaem** from CRAN. 
```{r}
install.packages("misaem")
 ```
### Using the misaem package
Basically,

1. `miss.glm` is the main function performing logistic regression with missing values.
2. `miss.lm` is the main function performing linear regression with missing values.

For more details, You can find the vignette, which illustrate the basic and further usage of misaem package:
```{r}
library(misaem)
vignette('misaem')
 ```

## Reference 
Logistic Regression with Missing Covariates
-- Parameter Estimation, Model Selection
and Prediction (2020, Jiang W., Josse J., Lavielle M., TraumaBase Group), [Computational Statistics & Data Analysis](https://doi.org/10.1016/j.csda.2019.106907).
