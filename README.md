# misaem
### Package R for "Stochastic Approximation EM for Logistic Regression with Missing Values (W. Jiang, J. Josse, M. Lavielle, 2018)"

``misaem`` is a method to apply statistical inference for logistic regression model with missing data. This method is based on likelihood, including 
1. A stochastic approximation version of EM algorithm based on Metropolis-Hasting sampling, to estimate the parameters of logistic regression;
2. Estimation of parameters' variance based one Louis formula;
3. Model selection procedure based on AIC or BIC.
