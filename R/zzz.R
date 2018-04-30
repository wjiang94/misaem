#' @importFrom utils packageDescription
#' @noRd
.onAttach <- function(libname, pkgname) {
  if (interactive()) {
    pdesc <- packageDescription(pkgname)
    packageStartupMessage('')
    packageStartupMessage(pdesc$Package,
                          " "
                          , pdesc$Version,
                          " par
                          "
                          ,pdesc$Author)
    packageStartupMessage(paste0('-> For help, type: help('
                                 ,pkgname,
                                 ')'))
    packageStartupMessage('')
  }}
