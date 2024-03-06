
.onAttach <- function(libname, pkgname) {
    packageStartupMessage("Attaching chiptsne2 version ",
                          utils::packageDescription("chiptsne2")$Version, ".")
}

