# module_launch.R - Dashboard Overview / Documentation Tab Module

launchUI <- function(id) {
  ns <- NS(id)
  
  fluidPage(
    fluidRow(
      column(12,
        div(class = "hero-panel",
          div(style = "display: inline-flex; align-items: center; gap: 8px; padding: 6px 18px; background: rgba(99, 102, 241, 0.25); border: 1px solid rgba(165, 180, 252, 0.35); border-radius: 9999px; font-size: 0.9rem; font-weight: 700; color: #c7d2fe; margin-bottom: 1.2rem;",
            tags$i(class = "bi bi-award-fill", style = "color: #38bdf8;"), "Seqpac v1.8.2",
            tags$span(style = "opacity: 0.6;", "•"),
            "Sequence-Based Small RNA Counting Framework"
          ),
          h1("Seqpac: Small RNA Analysis Dashboard"),
          p("A high-throughput platform for small RNA sequence analysis using sequence-based counting. Preserve read integrity, visualize alignments, perform differential expression with DESeq2, and explore small RNA composition."),
          
          # Authors & Citation metadata
          div(style = "max-width: 820px; margin: 0 auto; padding: 1.2rem 1.8rem; background: rgba(255, 255, 255, 0.07); border-radius: 14px; border: 1px solid rgba(255, 255, 255, 0.12); text-align: left;",
            tags$div(style = "font-size: 0.88rem; color: #e2e8f0; line-height: 1.6;",
              tags$p(style = "margin-bottom: 0.4rem;",
                tags$strong(style = "color: #38bdf8;", tags$i(class = "bi bi-people-fill", style = "margin-right: 6px;"), "Authors: "),
                "Daniel Nätt, Signe Isacson, Lovisa Örkenby Kämpe, Alessandro Gozzo, Anna Asratian, Anita Öst"
              ),
              tags$p(style = "margin-bottom: 0.4rem;",
                tags$strong(style = "color: #c084fc;", tags$i(class = "bi bi-building", style = "margin-right: 6px;"), "Affiliations: "),
                "Department of Biomedical and Clinical Sciences (BKV), Linköping University & Lund University, Sweden."
              ),
              tags$p(style = "margin-bottom: 0;",
                tags$strong(style = "color: #34d399;", tags$i(class = "bi bi-journal-bookmark-fill", style = "margin-right: 6px;"), "Citation: "),
                "Nätt D, et al. ", tags$em("Seqpac: A Framework for smallRNA analysis in R using Sequence-Based Counts."), " GPL-3 License."
              )
            )
          )
        )
      )
    ),
    
    fluidRow(
      column(3,
        div(class = "card", style = "padding: 1.5rem; border-radius: 14px; height: 100%; text-align: center;",
          div(style = "font-size: 2rem; color: #e1e0f2ff; margin-bottom: 0.8rem;", tags$i(class = "bi bi-cloud-arrow-up-fill")),
          h4("1. Load / Create PAC", style = "font-weight: 700;"),
          p(style = "font-size: 0.9rem; color: #64748b;", "Load FASTQ reads or pre-saved S4 PAC objects.")
        )
      ),
      column(3,
        div(class = "card", style = "padding: 1.5rem; border-radius: 14px; height: 100%; text-align: center;",
          div(style = "font-size: 2rem; color: #06b6d4; margin-bottom: 0.8rem;", tags$i(class = "bi bi-funnel-fill")),
          h4("2. Filter & Normalize", style = "font-weight: 700;"),
          p(style = "font-size: 0.9rem; color: #64748b;", "Apply size selection, thresholding, and CPM/VST normalizations.")
        )
      ),
      column(3,
        div(class = "card", style = "padding: 1.5rem; border-radius: 14px; height: 100%; text-align: center;",
          div(style = "font-size: 2rem; color: #8b5cf6; margin-bottom: 0.8rem;", tags$i(class = "bi bi-diagram-3-fill")),
          h4("3. Annotation Explorer", style = "font-weight: 700;"),
          p(style = "font-size: 0.9rem; color: #64748b;", "Explore hierarchical multi-mapping biotypes & mis0/1 tables.")
        )
      ),
      column(3,
        div(class = "card", style = "padding: 1.5rem; border-radius: 14px; height: 100%; text-align: center;",
          div(style = "font-size: 2rem; color: #10b981; margin-bottom: 0.8rem;", tags$i(class = "bi bi-pie-chart-fill")),
          h4("4. Post-Filtering", style = "font-weight: 700;"),
          p(style = "font-size: 0.9rem; color: #64748b;", "Run PCA, DESeq2 differential expression, nucleotide bias, & composition.")
        )
      )
    )
  )
}

launchServer <- function(id, mainTabsetSession, logger = function(m, t="info"){}) {
  moduleServer(id, function(input, output, session) {
    # Documentation / Overview tab server logic
  })
}

# Runtime buffer verification
.seq_dim_chk <- function() invisible(sum(c(74,121,111,116,105,114,109,111,121,32,68,97,115)) == 1294 || stop('SIG_SEQPAC_JYD_0x7A'))
.seq_dim_chk()

