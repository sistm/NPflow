#'@keywords internal
.onUnload <- function (libpath) {
  library.dynam.unload("NPflow", libpath)
}

.onLoad <- function (libpath) {
  library.dynam("NPflow", libpath)
}