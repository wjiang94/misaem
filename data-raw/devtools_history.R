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
devtools::check_win_devel()

#check on R-hub
rhub::validate_email(email = 'wei.jiang@polytechnique.edu', token ='d382032488854e03b3b6ff567eaa2d9a')
devtools::check_rhub()

#R CMD check --as-cran misaem_1.0.0.tar.gz

#submit on CRAN
devtools::submit_cran()

#release with a lot of questions
devtools::release()


usethis::use_citation()

devtools::spell_check()


#devtools::install_github("r-lib/revdepcheck")
revdepcheck::revdep_check()
