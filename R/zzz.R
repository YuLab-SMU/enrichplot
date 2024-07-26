##' @importFrom yulab.utils yulab_msg
.onAttach <- function(libname, pkgname) {
    options(check.tbl_tree.verbose = FALSE)
    packageStartupMessage(yulab_msg(pkgname))
}


