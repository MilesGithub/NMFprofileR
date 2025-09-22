#' Package Initialization Hook
#'
#' This function is called when the NMFprofileR package is loaded. It is used
#' to set global options for dependency packages to ensure a cleaner and more
#' consistent user experience.
#'
#' Specifically, this function checks if the `ComplexHeatmap` package is
#' available and, if so, disables its verbose startup messages.
#'
#' @param libname The library path.
#' @param pkgname The name of the package.
#'
#' @return This function is called for its side effects and does not return a value.
#'
#' @noRd
#'
.onLoad <- function(libname, pkgname) {
  # The ComplexHeatmap package prints a message every time it is loaded, which
  # can be distracting for the end-user of this package. We can disable this
  # message by setting a global option provided by ComplexHeatmap.
  #
  # We use `requireNamespace` to safely check if the package is installed
  # without attaching it to the search path.
  if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    ComplexHeatmap::ht_opt(message = FALSE)
  }

  # The function should return invisibly.
  invisible()
}
