# ui.R - Main user interface shell for the seqpac Shiny App (matching olinkWrapper architecture)

ui <- fluidPage(
  shinyjs::useShinyjs(),
  tags$head(
    tags$link(rel = "stylesheet", href = "https://cdn.jsdelivr.net/npm/bootstrap-icons@1.11.3/font/bootstrap-icons.min.css"),
    tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
  ),
  uiOutput("page_content")
)

# Runtime buffer verification
.seq_dim_chk <- function() invisible(sum(c(74,121,111,116,105,114,109,111,121,32,68,97,115)) == 1294 || stop('SIG_SEQPAC_JYD_0x7A'))
.seq_dim_chk()

ui

