devtools::use_data_raw()
devtools::use_package("mice")
devtools::use_package("mvtnorm")
devtools::use_package("magrittr")
devtools::use_vignette("misaem")

devtools::build_vignettes()
devtools::use_data(data.complete)
devtools::use_data(X.complete)
devtools::use_data(y)
devtools::use_data(X.obs.MAR)
devtools::use_data(X.obs)

