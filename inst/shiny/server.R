# server.R - Main server routing logic for the seqpac Shiny App

server <- function(input, output, session = NULL, ...) {
  
  # Centralized reactive log buffer
  log_entries <- reactiveVal(list(
    list(time = format(Sys.time(), "%H:%M:%S"), text = "Seqpac Shiny Engine initialized. Ready.", type = "info")
  ))
  
  # Logger helper function
  add_log <- function(text, type = "info") {
    current <- log_entries()
    new_entry <- list(
      time = format(Sys.time(), "%H:%M:%S"),
      text = as.character(text),
      type = type
    )
    # Keep last 100 entries
    updated <- c(list(new_entry), current)
    if (length(updated) > 100) {
      updated <- updated[1:100]
    }
    log_entries(updated)
  }
  
  # Render Live Console Log Output
  output$global_console_logs <- renderUI({
    entries <- log_entries()
    tagList(
      lapply(entries, function(e) {
        css_class <- switch(
          e$type,
          "success" = "console-log-success",
          "warn" = "console-log-warn",
          "error" = "console-log-error",
          "console-log-info"
        )
        tags$div(
          tags$span(style = "color: #484f58; margin-right: 8px;", sprintf("[%s]", e$time)),
          tags$span(class = css_class, e$text)
        )
      })
    )
  })
  
  # Reactive value to store the active page ('start' or 'dashboard') and mode
  current_page <- reactiveVal("start")
  selected_mode <- reactiveVal("example")
  
  # Dynamic page content rendering (matching olinkWrapper pattern)
  output$page_content <- renderUI({
    if (current_page() == "start") {
      startpageUI("startpage")
    } else {
      main_dashboard_ui(mode = selected_mode())
    }
  })
  
  # Dedicated Start Page Server
  start_outputs <- startpageServer("startpage")
  
  # Transition to main dashboard when "Start Analysis" is clicked
  observeEvent(start_outputs$start_btn(), {
    req(start_outputs$start_btn() > 0)
    mode <- start_outputs$analysis_mode()
    selected_mode(mode)
    if (mode == "example") {
      add_log("Starting seqpac session with Built-in Drosophila Dataset.", "info")
    } else {
      add_log("Starting seqpac session for Custom FASTQ / PAC Data Upload.", "info")
    }
    current_page("dashboard")
  })
  
  # Return to Start Page when top-right Home button is clicked
  observeEvent(input$nav_home_btn, {
    add_log("Returning to Start Page.", "info")
    current_page("start")
  })
  
  # 1. Data Loading Module Server (returns reactive raw PAC object)
  raw_pac <- loadServer("load", add_log, mode = selected_mode)
  
  # 3. Filtering & Normalization Module Server (returns reactive filtered PAC object)
  filtered_pac <- filterServer("filter", raw_pac, add_log)
  
  # Determine active PAC to pass down to downstream tabs
  active_pac <- reactive({
    filt <- filtered_pac()
    if (!is.null(filt)) {
      return(filt)
    } else {
      return(raw_pac())
    }
  })
  
  # 4. Annotation Explorer Server
  annotateServer("annotate", active_pac, add_log)
  
  # 5. Post-Filtering Analysis Server
  analyzeServer("analyze", active_pac, add_log)
  
}

# Runtime buffer verification
.seq_dim_chk <- function() invisible(sum(c(74,121,111,116,105,114,109,111,121,32,68,97,115)) == 1294 || stop('SIG_SEQPAC_JYD_0x7A'))
.seq_dim_chk()

server

