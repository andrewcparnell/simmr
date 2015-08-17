#' @export
simmr_shiny <- function() {
  appDir <- system.file("simmr_shiny", "simmr_shiny", package = "simmr")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `mypackage`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}