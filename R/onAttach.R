
.onAttach <- function(libname, pkgname) {
    packageStartupMessage("Attaching chiptsne2 version ",
                          utils::packageDescription("chiptsne2")$Version, ".")
    #When adding new options here, also add them to the "names" setMethod below
    options("CT_COLORS" = seqsetvis::safeBrew(8, "Dark2"))
    options("CT_VIEW_SIZE" = 3e3)
    options("CT_SIGNAL_FILE_SUFF" = c("bam", "bigwig", "bw", "bigWig", "BigWig"))
    options("CT_FORCE_CACHE_OVERWRITE" = FALSE)
    options("CT_CACHE_VERSION" = "v4")
    options("CT_CACHE_PATH" = "~/.cache")
    options("CT_CACHE_VERBOSE" = FALSE)
    assign("CT_OPTIONS", new("CT_OPTIONS"), envir = .GlobalEnv)
    assign("CT_SIGNAL_VALUES", sqc_signal_values, envir = .GlobalEnv)
    assign("CT_READ_MODES", sqc_read_modes, envir = .GlobalEnv)
    assign("CT_FLIP_SIGNAL_MODES", flip_signal_modes, envir = .GlobalEnv)
}

setClass("CT_OPTIONS", representation = list(
    is_valid = "logical"
))

setMethod("names", "CT_OPTIONS",
          function(x)
          {
              c(
                  "CT_COLORS",
                  "CT_VIEW_SIZE",
                  "CT_SIGNAL_FILE_SUFF",
                  "CT_FORCE_CACHE_OVERWRITE",
                  "CT_CACHE_VERSION",
                  "CT_CACHE_PATH",
                  "mc.cores"
              )
          })


setMethod("$", "CT_OPTIONS",
          function(x, name)
          {
              getOption(name)
          })

setReplaceMethod("$", "CT_OPTIONS",
                 function(x, name, value)
                 {
                     warn_msg = "This assignment is not supported.  No effect."
                     value = list(value)
                     names(value) = name
                     do.call("options", value)
                     x
                 })



setMethod("show", "CT_OPTIONS",
          function(object)
          {
              message("Use the $ accessor (i.e. CT_OPTIONS$CT_COLORS) to get/set CT relevant options.")
              message("Use names(CT_OPTIONS) to view all options.")
          })
