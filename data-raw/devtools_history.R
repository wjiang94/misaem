devtools::use_data_raw()
devtools::use_package("mice")
devtools::use_package("mvtnorm")

#vignette
devtools::use_package("magrittr")
devtools::use_vignette("misaem")
devtools::build_vignettes()

#check problems
devtools::check()

devtools::revdep_check()

#create a package bundle
devtools::build()

#check on windows and build a package bundle on windows
devtools::build_win()
#R CMD check --as-cran misaem_0.9.0.tar.gz


#submit on CRAN
devtools::submit_cran()
