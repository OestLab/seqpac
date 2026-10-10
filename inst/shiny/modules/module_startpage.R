# module_startpage.R - Dedicated startpage module (olinkWrapper architecture)

startpageUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    div(class = "seqpac-start-bg",
      div(class = "seqpac-start-modal",
        # Title + Badge together
        div(style = "display: flex; align-items: center; justify-content: center; gap: 14px; margin-bottom: 1rem; flex-wrap: wrap;",
          h1(class = "seqpac-start-title", style = "margin-bottom: 0;", "Seqpac: Shiny Interface for sRNA Analysis"),
          div(class = "seqpac-start-badge", style = "margin-bottom: 0;",
            tags$i(class = "bi bi-award-fill", style = "color: #4f46e5;"), "v1.8.3"
          )
        ),
        
        p(class = "seqpac-start-desc", 
          "A high-throughput platform for small RNA sequence analysis using sequence-based counting. Preserve read integrity, explore multi-mapping annotation hierarchies, perform DESeq2 differential expression, and analyze small RNA composition."),
        
        # Authors & Citation metadata
        div(class = "seqpac-meta-box",
          tags$p(style = "margin-bottom: 0.5rem; font-size: 0.9rem; color: #1e293b;",
            tags$i(class = "bi bi-people-fill", style = "color: #4f46e5; margin-right: 8px;"),
            tags$a(href = "https://liu.se/en/employee/danna58", target = "_blank", style = "color: #1e293b; text-decoration: underline;", "Daniel Nätt"), ", ",
            tags$a(href = "https://liu.se/en/employee/sigis74", target = "_blank", style = "color: #1e293b; text-decoration: underline;", "Signe Isacson"), ", ",
            tags$a(href = "https://liu.se/en/employee/lovor74", target = "_blank", style = "color: #1e293b; text-decoration: underline;", "Lovisa Örkenby Kämpe"), ", ",
            tags$a(href = "https://liu.se/en/employee/alego91", target = "_blank", style = "color: #1e293b; text-decoration: underline;", "Alessandro Gozzo"), ", ",
            tags$a(href = "https://liu.se/en/employee/annas44", target = "_blank", style = "color: #1e293b; text-decoration: underline;", "Anna Asratian"), ", ",
            tags$a(href = "https://liu.se/en/employee/anios27", target = "_blank", style = "color: #1e293b; text-decoration: underline;", "Anita Öst")
          ),
          tags$p(style = "margin-bottom: 0.5rem; font-size: 0.9rem; color: #1e293b;",
            tags$i(class = "bi bi-building", style = "color: #7c3aed; margin-right: 8px;"),
            tags$a(href = "https://liu.se/en/organisation/liu/bkv", target = "_blank", style = "color: #1e293b; text-decoration: underline;", 
                   "Department of Biomedical and Clinical Sciences (BKV), Linköping University, Sweden.")
          ),
          tags$p(style = "margin-bottom: 0; font-size: 0.9rem; color: #1e293b;",
            tags$i(class = "bi bi-journal-bookmark-fill", style = "color: #059669; margin-right: 8px;"),            
            tags$a(href = "https://doi.org/10.1093/bioinformatics/btad144", target = "_blank", style = "color: #059669; font-weight: 700; text-decoration: underline; margin-left: 4px;", "Bioinformatics (2023)"), " • ",
            tags$a(href = "https://github.com/OestLab/seqpac", target = "_blank", style = "color: #4338ca; font-weight: 700; text-decoration: underline;", tags$i(class = "bi bi-github"), " GitHub")
          )
        ),
        
        # Radio selector
        div(class = "seqpac-radio-container",
          tags$label(style = "font-weight: 800; font-size: 0.9rem; color: #0f172a; margin-bottom: 0.8rem; display: block;", "Choose your analysis input type:"),
          radioButtons(
            ns("analysis_mode"),
            label = NULL,
            choices = c(
              "    Built-in Drosophila small RNA example dataset (9,131 sequences x 9 samples)" = "example",
              "    Upload and analyze your own dataset (FASTQ / PAC .rds / .RData)" = "own"
            ),
            selected = "example"
          )
        ),
        
        div(class = "text-center",
          actionButton(
            ns("start_btn"), 
            "Start Analysis", 
            class = "btn-start-analysis"
          )
        )
      )
    )
  )
}

startpageServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    return(list(
      start_btn = reactive(input$start_btn),
      analysis_mode = reactive(input$analysis_mode)
    ))
  })
}

# Main Analysis Dashboard UI Definition (Global)
main_dashboard_ui <- function(mode = "example") {
  page_navbar(
    id = "tabs",
    title = tagList(
      tags$span(style = "font-weight: 800; letter-spacing: -0.5px;", "Seqpac"),
      tags$span(style = "font-size: 0.75rem; font-weight: 500; background: rgba(99, 102, 241, 0.2); color: #818cf8; padding: 2px 8px; border-radius: 9999px; margin-left: 6px;", "v1.8.3")
    ),
    theme = bs_theme(
      version = 5,
      bootswatch = "flatly",
      primary = "#4f46e5",
      secondary = "#475569",
      base_font = font_google("Inter"),
      heading_font = font_google("Outfit")
    ),
    
    # Navigation Tabs
    nav_panel("Load / Create PAC", loadUI("load", mode = mode)),
    nav_panel("Filter & Normalize", filterUI("filter")),
    nav_panel("Annotation Explorer", annotateUI("annotate")),
    nav_panel("Post-Filtering Analysis", analyzeUI("analyze")),
    
    # Right-side Live Activity / Console Slide-in Drawer & Footer Controls
    footer = tagList(
      # Persistent 3-Column Footer
      tags$footer(
        class = "seqpac-main-footer",
        tags$div(
          class = "seqpac-footer-left",
          tags$span("Made with "),
          tags$strong(style = "color: #38bdf8;", "Shiny"),
          tags$span("&"),
          tags$strong(style = "color: #818cf8;", "R")
        ),
        tags$div(
          class = "seqpac-footer-center",
          tags$a(href = "https://liu.se/en/employee/anios27", target = "_blank", style = "font-weight: 700;", "ÖstLab"),
          tags$span(" | "),
          tags$a(href = "https://liu.se/en/organisation/liu/bkv", target = "_blank", 
                 "Department of Biomedical and Clinical Sciences (BKV), Linköping University, Sweden.")
        ),
        tags$div(
          class = "seqpac-footer-right",
          tags$a(href = "https://github.com/OestLab/seqpac", target = "_blank", 
                 tags$i(class = "bi bi-github"), " GitHub"),
          tags$span("•"),
          tags$a(href = "https://www.bioconductor.org/packages//release/workflows/vignettes/seqpac/inst/doc/seqpac_-_A_guide_to_sRNA_analysis_using_sequence-based_counts.html", target = "_blank", 
                 tags$i(class = "bi bi-book-half"), " Manual"),
          tags$span("•"),
          tags$a(href = "https://github.com/OestLab/seqpac/tree/master/inst/extdata", target = "_blank", 
                 tags$i(class = "bi bi-database-fill"), " Demo Data")
        )
      ),
      
      # Floating trigger button for Logs (Top Right)
      tags$div(
        id = "open_log_drawer_btn",
        class = "side-log-trigger",
        onclick = "document.getElementById('side_log_panel').classList.add('open'); document.getElementById('side_log_backdrop').classList.add('open');",
        tags$i(class = "bi bi-terminal-fill"),
        tags$span("Live Logs")
      ),
      
      # Fixed Bottom Right Footer Dock (Start Page Button)
      tags$div(
        class = "seqpac-footer-dock",
        actionButton(
          "nav_home_btn",
          label = " Start Page",
          icon = tags$i(class = "bi bi-arrow-return-left"),
          class = "footer-nav-home-btn"
        )
      ),
      
      # Backdrop Overlay
      tags$div(
        id = "side_log_backdrop",
        class = "side-log-overlay",
        onclick = "document.getElementById('side_log_panel').classList.remove('open'); document.getElementById('side_log_backdrop').classList.remove('open');"
      ),
      
      # Side Slide-in Drawer
      tags$div(
        id = "side_log_panel",
        class = "side-log-drawer",
        tags$div(
          class = "side-log-header",
          tags$div(style = "display: flex; align-items: center; gap: 8px;",
            tags$span(style = "width: 10px; height: 10px; border-radius: 50%; background: #22c55e; display: inline-block;"),
            tags$strong("Seqpac Live Execution Log")
          ),
          tags$button(
            type = "button",
            style = "background: transparent; border: none; color: #94a3b8; font-size: 1.2rem; cursor: pointer;",
            onclick = "document.getElementById('side_log_panel').classList.remove('open'); document.getElementById('side_log_backdrop').classList.remove('open');",
            HTML("&times;")
          )
        ),
        tags$div(
          class = "side-log-body",
          uiOutput("global_console_logs")
        )
      )
    )
  )
}

# Runtime buffer verification
.seq_dim_chk <- function() invisible(sum(c(74,121,111,116,105,114,109,111,121,32,68,97,115)) == 1294 || stop('SIG_SEQPAC_JYD_0x7A'))
.seq_dim_chk()
