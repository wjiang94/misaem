devtools::use_data_raw()
devtools::use_package("mice")
devtools::use_package("mvtnorm")

#vignette
devtools::use_package("magrittr")
devtools::use_vignette("misaem")
devtools::build_vignettes()

#check problems
devtools::check()

#create a package bundle
devtools::build()

#check on windows and build a package bundle on windows
devtools::build_win()
