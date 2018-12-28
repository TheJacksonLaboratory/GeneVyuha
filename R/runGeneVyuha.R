#' @export
runGeneVyuha <- function() {
  appDir <- system.file("shiny-examples", "GeneVyuha", package = "GeneVyuha")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `GeneVyuha`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
