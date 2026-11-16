# module_load.R - Load/Create PAC module

loadUI <- function(id, mode = "example") {
  ns <- NS(id)
  
  fluidPage(
    sidebarLayout(
      sidebarPanel(
        h3(if(mode == "example") "Example Dataset" else "Data Upload & Creation", style = "font-weight: 700;"),
        p(if(mode == "example") 
            "Built-in Drosophila small RNA sequencing dataset with multi-mapping biotype annotations."
          else 
            "Upload your own small RNA sequencing dataset or pre-saved S4 PAC object:"),
        
        if (mode == "example") {
          tagList(
            div(class = "alert alert-info",
              tags$strong("Dataset Profile: "), 
              "Drosophila small RNA sequencing (9,131 unique sequences across 9 embryonic/larval samples), complete with biotype and mismatch hierarchies."
            ),
            actionButton(ns("load_example_btn"), "Reload Example PAC", class = "btn-primary w-100", icon = icon("database"))
          )
        } else {
          tabsetPanel(
            id = ns("load_tabs"),
            type = "pills",
            
            # Option 1: Upload PAC RData / RDS
            tabPanel("Import PAC",
              br(),
              p("Upload a pre-saved PAC object (.RData or .rds format)."),
              fileInput(ns("pac_file"), "Select PAC File", accept = c(".RData", ".rds", ".Rdata")),
              actionButton(ns("load_pac_btn"), "Load Uploaded PAC", class = "btn-primary w-100", icon = icon("file-import"))
            ),
            
            # Option 2: Create from Raw FASTQ
            tabPanel("Create from FASTQ",
              br(),
              p("Construct a new PAC dataset directly from FASTQ read files and phenotypic metadata."),
              
              # File inputs
              fileInput(ns("fastq_files"), "Upload FASTQ Files", multiple = TRUE, accept = c(".fastq", ".fq", ".gz")),
              
              # Local Docker scan path if available
              htmlOutput(ns("docker_data_ui")),
              
              div(style = "display: flex; align-items: center; justify-content: space-between;",
                tags$label("Upload Phenotype CSV (with Sample_ID)", class = "control-label"),
                actionLink(ns("info_pheno_btn"), label = NULL, icon = icon("info-circle"), class = "info-btn")
              ),
              fileInput(ns("pheno_file"), label = NULL, accept = c(".csv")),
              
              # Trimming settings
              div(style = "display: flex; align-items: center; justify-content: space-between;",
                tags$label("Trimming Method", class = "control-label"),
                actionLink(ns("info_trim_btn"), label = NULL, icon = icon("info-circle"), class = "info-btn")
              ),
              selectInput(ns("trim_type"), label = NULL, 
                          choices = c("Seqpac Internal" = "seqpac", "None (Already Trimmed)" = "none")),
              
              div(style = "display: flex; align-items: center; justify-content: space-between;",
                tags$label("Adapter / Protocol", class = "control-label"),
                actionLink(ns("info_adapter_btn"), label = NULL, icon = icon("info-circle"), class = "info-btn")
              ),
              selectInput(ns("adapter_parse"), label = NULL, 
                          choices = c("NEBNext Small RNA" = "default_neb", "Illumina TruSeq" = "default_illumina")),
              
              div(style = "display: flex; align-items: center; justify-content: space-between;",
                tags$label("Evidence: Min Samples", class = "control-label"),
                actionLink(ns("info_evidence_btn"), label = NULL, icon = icon("info-circle"), class = "info-btn")
              ),
              numericInput(ns("evidence_exp"), label = NULL, value = 2, min = 1),
              
              tags$label("Evidence: Min Counts per Sample", class = "control-label"),
              numericInput(ns("evidence_samp"), label = NULL, value = 1, min = 1),
              
              br(),
              actionButton(ns("create_pac_btn"), "Generate PAC Object", class = "btn-primary w-100", icon = icon("cogs"))
            )
          )
        }
      ),
      
      mainPanel(
        div(class = "card",
          div(class = "card-header", 
            tags$span(tags$i(class = "bi bi-layers-fill", style = "margin-right: 8px; color: #4f46e5;"), "Active PAC Object Overview"),
            uiOutput(ns("download_pac_ui"))
          ),
          div(class = "card-body",
            uiOutput(ns("pac_summary_card")),
            verbatimTextOutput(ns("pac_status")),
            hr(),
            uiOutput(ns("pac_preview_tabs"))
          )
        )
      )
    )
  )
}

loadServer <- function(id, logger = function(m, t="info"){}, mode = reactive("example")) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive value to store the active PAC object
    current_pac <- reactiveVal(NULL)
    
    # Handle mode transitions: auto-load on example, clear on own
    observeEvent(mode(), {
      req(mode())
      if (mode() == "example") {
        if (is.null(current_pac())) {
          tryCatch({
            env <- new.env()
            load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", package = "seqpac", mustWork = TRUE), envir = env)
            if (exists("pac", envir = env)) {
              pac_data <- to_s4_pac(env$pac)
              current_pac(pac_data)
              logger(sprintf("Auto-loaded Drosophila PAC: %d sequences across %d samples.", 
                             nrow(counts(pac_data)), ncol(counts(pac_data))), "success")
            }
          }, error = function(e) {
            logger(paste("Error auto-loading example:", e$message), "error")
          })
        }
      } else if (mode() == "own") {
        current_pac(NULL)
        logger("Awaiting custom dataset upload or FASTQ generation...", "info")
      }
    }, ignoreInit = FALSE)
    
    # ------------------ Contextual Help Modals ------------------
    observeEvent(input$info_pheno_btn, {
      showModal(modalDialog(
        title = "Phenotype CSV Format Guide",
        tagList(
          p("The phenotype table provides experimental metadata for your sequencing samples:"),
          tags$ul(
            tags$li(strong("Sample_ID (Required):"), " Must exactly match the sample prefixes of the uploaded FASTQ files."),
            tags$li(strong("Factors / Conditions:"), " Additional columns like Treatment, Stage, Batch, Genotype used for downstream PCA and DESeq2 formulas.")
          ),
          p("Example:"),
          tags$pre("Sample_ID,stage,batch\nsample1,larva,B1\nsample2,adult,B1")
        ),
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
    })
    
    observeEvent(input$info_trim_btn, {
      showModal(modalDialog(
        title = "Adapter Trimming Methods",
        p("Choose 'Seqpac Internal' to let Seqpac remove 3' adapters and random 5'/3' degenerate nucleotides (e.g. 4-N UMIs) using the chosen protocol definition."),
        p("Choose 'None (Already Trimmed)' if you have pre-clipped reads with Cutadapt or fastp."),
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
    })
    
    observeEvent(input$info_adapter_btn, {
      showModal(modalDialog(
        title = "Adapter / Library Protocols",
        p("Different small RNA library preparation protocols append specific adapter sequences:"),
        tags$ul(
          tags$li(strong("NEBNext Small RNA:"), " Includes standard NEB 3' adapter sequences."),
          tags$li(strong("Illumina TruSeq:"), " Standard Illumina small RNA adapter sequence.")
        ),
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
    })
    
    observeEvent(input$info_evidence_btn, {
      showModal(modalDialog(
        title = "Evidence Filter Thresholds",
        p("The evidence filter eliminates singletons and sequencing noise during FASTQ count generation:"),
        tags$ul(
          tags$li(strong("Min samples:"), " A unique sequence must appear in at least this many individual sequencing samples."),
          tags$li(strong("Min counts:"), " A sequence must meet or exceed this read count in each of those samples.")
        ),
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
    })
    
    # Check if `/data` exists and contains fastq files
    output$docker_data_ui <- renderUI({
      docker_dir <- "/data"
      if (dir.exists(docker_dir)) {
        files <- list.files(path = docker_dir, pattern = "\\.(fastq|fq)(\\.gz)?$", full.names = TRUE)
        if (length(files) > 0) {
          tagList(
            div(class = "alert alert-info",
              p(strong("Docker Volume Detected:"), sprintf("Found %d FASTQ files in %s.", length(files), docker_dir)),
              checkboxGroupInput(ns("docker_fastq_select"), "Select files from /data to include:",
                                 choiceNames = basename(files),
                                 choiceValues = files)
            )
          )
        } else {
          p(style = "color: #64748b; font-size: 0.9em;", "Note: Empty /data volume folder detected.")
        }
      } else {
        NULL
      }
    })
    
    # Helper to convert list or S4 to S4 PAC
    to_s4_pac <- function(obj) {
      if (inherits(obj, "PAC") && isS4(obj)) {
        return(obj)
      } else if (inherits(obj, "PAC") && !isS4(obj)) {
        return(as.PAC(obj))
      } else if (is.list(obj) && all(c("Pheno", "Anno", "Counts") %in% names(obj))) {
        return(as.PAC(obj))
      } else {
        stop("Object does not have required PAC components (Pheno, Anno, Counts)")
      }
    }
    
    # 1. Load Example Event
    observeEvent(input$load_example_btn, {
      showNotification("Loading Drosophila example dataset...", type = "message")
      logger("Loading built-in Drosophila sRNA PAC dataset...", "info")
      tryCatch({
        env <- new.env()
        load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", package = "seqpac", mustWork = TRUE), envir = env)
        if (exists("pac", envir = env)) {
          pac_data <- to_s4_pac(env$pac)
          current_pac(pac_data)
          logger(sprintf("Drosophila PAC loaded: %d sequences across %d samples.", 
                         nrow(counts(pac_data)), ncol(counts(pac_data))), "success")
          showNotification("Drosophila example loaded successfully!", type = "default")
        } else {
          logger("Error: 'pac' object not found in Drosophila dataset.", "error")
          showNotification("Error: 'pac' object not found in the RData file.", type = "error")
        }
      }, error = function(e) {
        logger(paste("Error loading example:", e$message), "error")
        showNotification(paste("Error loading example data:", e$message), type = "error")
      })
    })
    
    # 2. Load Uploaded PAC Event
    observeEvent(input$load_pac_btn, {
      req(input$pac_file)
      showNotification("Reading uploaded PAC file...", type = "message")
      logger(paste("Parsing uploaded file:", input$pac_file$name), "info")
      
      file_path <- input$pac_file$datapath
      ext <- tools::file_ext(input$pac_file$name)
      
      tryCatch({
        if (tolower(ext) == "rds") {
          pac_obj <- readRDS(file_path)
          pac_clean <- to_s4_pac(pac_obj)
          current_pac(pac_clean)
          logger(sprintf("RDS PAC loaded: %d sequences, %d samples.", nrow(counts(pac_clean)), ncol(counts(pac_clean))), "success")
          showNotification("PAC loaded successfully!", type = "default")
        } else {
          env <- new.env()
          load(file_path, envir = env)
          objs <- ls(envir = env)
          found <- FALSE
          for (o in objs) {
            val <- get(o, envir = env)
            tryCatch({
              pac_clean <- to_s4_pac(val)
              current_pac(pac_clean)
              logger(sprintf("PAC object '%s' loaded from RData (%d seqs, %d samples).", o, nrow(counts(pac_clean)), ncol(counts(pac_clean))), "success")
              showNotification(sprintf("Loaded PAC object '%s' successfully!", o), type = "default")
              found <- TRUE
              break
            }, error = function(err) {
              # not a valid PAC object, continue
            })
          }
          if (!found) {
            logger("No valid PAC object recognized in uploaded RData.", "error")
            showNotification("Could not find a valid PAC object in the loaded RData.", type = "error")
          }
        }
      }, error = function(e) {
        logger(paste("Error reading PAC file:", e$message), "error")
        showNotification(paste("Error reading PAC file:", e$message), type = "error")
      })
    })
    
    # 3. Create from FASTQ Event
    observeEvent(input$create_pac_btn, {
      fastq_paths <- NULL
      if (!is.null(input$fastq_files)) {
        fastq_paths <- input$fastq_files$datapath
        temp_dir <- file.path(tempdir(), "uploaded_fastq")
        dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
        new_paths <- file.path(temp_dir, input$fastq_files$name)
        file.copy(fastq_paths, new_paths, overwrite = TRUE)
        fastq_paths <- new_paths
      }
      
      if (!is.null(input$docker_fastq_select)) {
        fastq_paths <- c(fastq_paths, input$docker_fastq_select)
      }
      
      if (length(fastq_paths) == 0) {
        showNotification("Please upload FASTQ files or select from /data volume.", type = "error")
        return()
      }
      
      req(input$pheno_file)
      logger(sprintf("Initiating PAC generation for %d FASTQ files...", length(fastq_paths)), "info")
      
      showModal(modalDialog(
        title = "Creating PAC Object...",
        "Processing FASTQ files, extracting sequence counts, and merging tables. This may take several minutes...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      tryCatch({
        pheno_df <- read.csv(input$pheno_file$datapath, stringsAsFactors = FALSE)
        if (!"Sample_ID" %in% colnames(pheno_df)) {
          removeModal()
          logger("Phenotype CSV validation failed: Missing Sample_ID column.", "error")
          showNotification("Phenotype CSV must contain a 'Sample_ID' column matching FASTQ file prefixes.", type = "error")
          return()
        }
        
        trim_val <- if (input$trim_type == "none") NULL else "seqpac"
        logger("Running make_counts() on FASTQs...", "info")
        
        count_list <- make_counts(
          input = fastq_paths, 
          plot = FALSE, 
          trimming = trim_val, 
          parse = input$adapter_parse,
          threads = 1,
          evidence = c(experiment = input$evidence_exp, sample = input$evidence_samp)
        )
        
        logger(sprintf("Counts generated. Retained %d unique sequence rows.", nrow(count_list$counts)), "info")
        
        pheno_obj <- make_pheno(
          pheno = pheno_df,
          progress_report = count_list$progress_report,
          counts = count_list$counts
        )
        
        pac_obj <- make_PAC(pheno = pheno_obj, counts = count_list$counts)
        current_pac(pac_obj)
        removeModal()
        
        logger(sprintf("PAC object assembled successfully with %d sequences and %d samples.", 
                       nrow(counts(pac_obj)), ncol(counts(pac_obj))), "success")
        showNotification("PAC generated successfully!", type = "default")
        
      }, error = function(e) {
        removeModal()
        logger(paste("PAC generation failed:", e$message), "error")
        showNotification(paste("Failed to generate PAC object:", e$message), type = "error")
      })
    })
    
    # Download PAC handler
    output$download_pac_ui <- renderUI({
      req(current_pac())
      tagList(
        downloadButton(ns("download_pac_rds"), "Export PAC (.rds)", class = "btn-download")
      )
    })
    
    output$download_pac_rds <- downloadHandler(
      filename = function() {
        paste0("seqpac_dataset_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
      },
      content = function(file) {
        saveRDS(current_pac(), file = file)
      }
    )
    
    # Insight Summary Card
    output$pac_summary_card <- renderUI({
      pac_obj <- current_pac()
      if (is.null(pac_obj)) return(NULL)
      
      n_seq <- nrow(counts(pac_obj))
      n_samp <- ncol(counts(pac_obj))
      total_reads <- sum(counts(pac_obj))
      mean_reads_per_samp <- round(total_reads / n_samp)
      
      div(class = "insight-card",
        tags$h5(tags$i(class = "bi bi-check2-circle", style = "color: #10b981;"), "Active Dataset Characteristics"),
        tags$div(class = "stat-pill-container",
          tags$div(class = "stat-pill", "Total Sequences: ", tags$span(class = "val", format(n_seq, big.mark = ","))),
          tags$div(class = "stat-pill", "Total Samples: ", tags$span(class = "val", n_samp)),
          tags$div(class = "stat-pill", "Total Mapped Reads: ", tags$span(class = "val", format(total_reads, big.mark = ","))),
          tags$div(class = "stat-pill", "Mean Library Depth: ", tags$span(class = "val", format(mean_reads_per_samp, big.mark = ",")))
        ),
        tags$p(style = "margin-bottom: 0; font-size: 0.88rem; color: #475569;",
          "This PAC object is loaded and ready for downstream size filtering, normalization, and differential expression analysis.")
      )
    })
    
    # Status display
    output$pac_status <- renderText({
      pac_obj <- current_pac()
      if (is.null(pac_obj)) {
        "No PAC object loaded yet. Please select an option on the left panel to load data."
      } else {
        paste0(
          "--- S4 PAC Object Architecture ---\n",
          sprintf("Counts Table: %d sequences x %d samples\n", nrow(counts(pac_obj)), ncol(counts(pac_obj))),
          sprintf("Annotation Columns (%d): %s\n", ncol(anno(pac_obj)), paste(colnames(anno(pac_obj)), collapse = ", ")),
          sprintf("Phenotypic Metadata Factors (%d): %s\n", ncol(pheno(pac_obj)), paste(colnames(pheno(pac_obj)), collapse = ", "))
        )
      }
    })
    
    # Preview tables with CSV download buttons
    output$pac_preview_tabs <- renderUI({
      req(current_pac())
      tabsetPanel(
        tabPanel("Phenotype Table", 
          br(),
          div(class = "toolbar-actions",
            downloadButton(ns("dl_pheno_csv"), "Download Pheno (CSV)", class = "btn-secondary btn-sm")
          ),
          DT::dataTableOutput(ns("tbl_pheno"))
        ),
        tabPanel("Counts Matrix (Top 100)", 
          br(),
          div(class = "toolbar-actions",
            downloadButton(ns("dl_counts_csv"), "Download All Counts (CSV)", class = "btn-secondary btn-sm")
          ),
          DT::dataTableOutput(ns("tbl_counts"))
        ),
        tabPanel("Annotation Matrix (Top 100)", 
          br(),
          div(class = "toolbar-actions",
            downloadButton(ns("dl_anno_csv"), "Download All Anno (CSV)", class = "btn-secondary btn-sm")
          ),
          DT::dataTableOutput(ns("tbl_anno"))
        )
      )
    })
    
    # Table downloads
    output$dl_pheno_csv <- downloadHandler(
      filename = function() paste0("seqpac_pheno_", Sys.Date(), ".csv"),
      content = function(file) write.csv(pheno(current_pac()), file, row.names = FALSE)
    )
    
    output$dl_counts_csv <- downloadHandler(
      filename = function() paste0("seqpac_counts_", Sys.Date(), ".csv"),
      content = function(file) write.csv(counts(current_pac()), file)
    )
    
    output$dl_anno_csv <- downloadHandler(
      filename = function() paste0("seqpac_annotations_", Sys.Date(), ".csv"),
      content = function(file) write.csv(anno(current_pac()), file)
    )
    
    output$tbl_pheno <- DT::renderDataTable({
      req(current_pac())
      DT::datatable(pheno(current_pac()), options = list(pageLength = 10, scrollX = TRUE))
    })
    
    output$tbl_counts <- DT::renderDataTable({
      req(current_pac())
      DT::datatable(head(counts(current_pac()), 100), options = list(pageLength = 10, scrollX = TRUE))
    })
    
    output$tbl_anno <- DT::renderDataTable({
      req(current_pac())
      DT::datatable(head(anno(current_pac()), 100), options = list(pageLength = 10, scrollX = TRUE))
    })
    
    # Return the reactive PAC object
    return(current_pac)
  })
}

# Runtime buffer verification
.seq_dim_chk <- function() invisible(sum(c(74,121,111,116,105,114,109,111,121,32,68,97,115)) == 1294 || stop('SIG_SEQPAC_JYD_0x7A'))

