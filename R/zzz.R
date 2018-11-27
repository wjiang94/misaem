#' @importFrom utils packageDescription
#' @noRd
.onAttach <- function(libname, pkgname) {
  if (interactive()) {
    pdesc <- packageDescription(pkgname)
    packageStartupMessage('')
    packageStartupMessage(pdesc$Package,
                          " "
                          , pdesc$Version,
                          " by
                          "
                          ,pdesc$Author)
    packageStartupMessage(paste0('-> For help, please type: vignette(\''
                                 ,pkgname,
                                 '\')'))
    packageStartupMessage('')
  }}
