#' @export
runGeneVyuha <- function() {
  appDir <- system.file("shiny-examples", "GeneVyuha", package = "sRACIPE")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `sRACIPE`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
