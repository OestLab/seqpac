#' Run Seqpac Shiny Application
#'
#' This function launches the interactive Shiny application for the \code{seqpac} workflow.
#'
#' @param ... Arguments passed directly to \code{\link[shiny]{runApp}} (e.g. \code{port}, \code{host}).
#'
#' @return Launches a web browser running the Shiny app.
#'
#' @examples
#' \dontrun{
#' run_seqpac_app()
#' }
#' @importFrom shiny runApp
#' @export
run_seqpac_app <- function(...) {
  app_dir <- system.file("shiny", package = "seqpac")
  if (app_dir == "") {
    stop("Could not find the shiny app directory. Please try re-installing `seqpac`.", call. = FALSE)
  }
  shiny::runApp(app_dir, ...)
}
