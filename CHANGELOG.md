# Changelog

All notable changes to the `seqpac` package will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [0.1.0] – 2026-08-18

### Added
- **Dedicated Start Page Architecture (`module_startpage.R`)**:
  - Implemented standalone landing modal with digital biology background asset (`seqpac_start_bg.jpg`).
  - Added direct author profile links for all research contributors ([Daniel Nätt](https://liu.se/en/employee/danna58), [Signe Isacson](https://liu.se/en/employee/sigis74), [Lovisa Örkenby Kämpe](https://liu.se/en/employee/lovor74), [Alessandro Gozzo](https://liu.se/en/employee/alego91), [Anna Asratian](https://liu.se/en/employee/annas44), and [Anita Öst](https://liu.se/en/employee/anios27)).
  - Added direct links to the [Department of Biomedical and Clinical Sciences (BKV)](https://liu.se/en/organisation/liu/bkv), [Bioinformatics Paper (DOI)](https://doi.org/10.1093/bioinformatics/btad144), and [GitHub Repository](https://github.com/OestLab/seqpac).
  - Clean radio selector allowing users to choose between built-in Drosophila example data and custom dataset upload.
- **Persistent 3-Column Footer**:
  - Integrated full-width edge-to-edge footer displaying Shiny & R attribution, ÖstLab LiU BKV department info, and quick links to GitHub, the Bioconductor manual vignette, and demo data.
- **Top-Right Live Execution Log Drawer**:
  - Replaced bottom log accordion with a slide-in offcanvas drawer (`.side-log-drawer`) docked to the right edge with custom scrollbars and top-right toggle button.
- **Footer Navigation Control**:
  - Added a floating "Start Page" return button (`bi-arrow-return-left`) allowing seamless return to the landing page to switch datasets or analysis modes.

### Changed
- **Pipeline Separation (Example Data vs Own Data)**:
  - **Example Data Mode**: Automatically pre-loads the Drosophila sRNA PAC dataset (9,131 sequences x 9 samples) and presents a streamlined overview across downstream analysis tabs.
  - **Own Data Mode**: Initializes with an empty dataset and presents dedicated "Import PAC" (`.rds` / `.RData`) and "Create from FASTQ" pipelines.
- **Solid High-Contrast Theme & Standard Web Proportions (`style.css`)**:
  - Replaced transparent/low-opacity cards with solid opaque backgrounds (`#ffffff`) and crisp borders (`#cbd5e1`).
  - Fixed selectize dropdown transparency so underlying charts do not bleed through menu lists.
  - Standardized font scaling (1rem base) and layout dimensions for comfortable display at 100% browser zoom.
  - Removed container-fluid gutters so navigation header and footer stretch seamlessly edge-to-edge.
  - Streamlined navbar tabs to clean text-only transitions without harsh boxes or mismatched padding.

### Fixed
- **S4 PAC `.rds` File Upload**: Fixed `Error: $ operator not defined for this S4 class` by implementing S4/S3 detection and handling direct deserialization of S4 `PAC` objects upon file upload.
- **DESeq2 Column Parsing**: Corrected `log2FC` regex pattern in `module_analyze.R` to ensure up/down-regulated sequence counts calculate accurately.

---

## [0.0.1] – 2026-06-24

### Added
- **Shiny Application**: Integrated an interactive, modular Shiny application under `inst/shiny/` to run the small RNA sequence counts analysis workflow described in the vignette:
  - **Launch Page** (`module_launch.R`): A landing page introducing the workflow steps.
  - **Load / Create PAC** (`module_load.R`): Supports loading the Drosophila example PAC, uploading custom `.RData`/`.rds` files, or generating a PAC from raw FASTQ files and a pheno CSV.
  - **Filter & Normalize** (`module_filter.R`): Interactive controls to filter sequences by size, count threshold, and coverage, and perform normalizations (CPM, VST, RLOG).
  - **Annotation Explorer** (`module_annotate.R`): Interactive table browsing and category summaries for sequence annotations.
  - **Post-Filtering Analysis** (`module_analyze.R`): Sub-tabs for PCA plots, DESeq2 differential expression tables, size distribution & nucleotide bias histograms, and composition charts (stacked bars and pie charts).
  - **Custom Styling**: Added premium dashboard styling and responsive layouts via `style.css` (Inter and Outfit google fonts).
- **Package Launcher**: Added `run_seqpac_app()` in `R/run_app.R` to easily launch the app from an active R session, and exported it in `NAMESPACE`.
- **Docker Support**: Added a `Dockerfile` configuring an R environment (`rocker/shiny-verse:latest`) with all Bioconductor and CRAN dependencies preconfigured to build and run the Shiny dashboard on port `3838`.
- **Instructions**: Added `run-shiny.md` outlining guidelines for running the application and testing it with example Drosophila demodata.

### Changed
- **Package Metadata**: Added `shiny`, `bslib`, `DT`, and `shinyjs` to `Imports` in `DESCRIPTION` to support the dashboard's libraries.

---

## [1.8.2] – 2026-04-30

### Changed
- **README**: Updated installation instructions to replace the deprecated `devtools::install_github()` call with the recommended `BiocManager::install()` approach.

---

## [1.8.1] – 2026-03-31

### Fixed
- **ggplot2 compatibility** (PR #54): Resolved deprecation warnings introduced by ggplot2 ≥ 3.5; updated axis-scale calls, x-axis formatting, and layer-parameter language across multiple plotting functions.
- **Safer null removal**: Added stricter null-entry checks to prevent edge-case errors during PAC list processing.
- **`tRNA_class.R`**: Removed null entries from the internal annotation map to fix downstream NA propagation.
- **`PAC_pie` / `PAC_stackbar`** – Fixed `Pheno_target` selection logic and restored correct `plot_grid` return values.
- **Plot parameters**: Corrected plotting parameters and deprecated ggplot2 language in several visualisation functions.
- **User-facing messages**: Standardised "Script terminated by user." wording to be consistent across functions.

---

## [1.8.0] – 2026-03-25

### Added
- **Documentation & examples** (PR #53): Expanded vignette and help pages for `make_conv` and `PAC_analyze` with worked examples and updated explanatory text.
- **`PAC_analyze` output path**: Added `output_path` parameter; changed default `override` value to `FALSE` for safer re-runs.

### Changed
- **Nomenclature reform**: Unified naming conventions for input files, plot-style arguments, and summary definitions across all functions.
  - `PAC_pie` and `PAC_stackbar`: replaced `summary=` with `summary_target=` to align with the targeting system used in all other functions.
- **`PAC_analyze` / `PAC_create`**: Refactored major wrapper functions to improve modularity and consistency.
- **Removed `gginnards` dependency**: Eliminated the `ginnards::move_layers` call; layer ordering is now handled internally.

### Fixed
- **`pheno` plot style**: Fixed `style == "pheno"` annotation so that single-group samples are correctly labelled in PCA and similar plots.

---

## [1.6.0] – 2025-09-26

### Added
- **Annotation workflow quick guide** (PR #47): Added a dedicated quick-start guide to the vignette covering the end-to-end annotation workflow.
- **`make_trim` file-extension support**: Extended `make_trim` to recognise `.fq` in addition to `.fastq` file endings, preventing silent failures on non-standard naming.

### Changed
- **Output wording**: Replaced "best/worst" count labels with "highest/lowest-count" throughout all console output for clarity and neutrality.
- **ggplot2 modernisation**: Updated `fviz_pca` parsing and axis-scale calls; replaced base-R plots with ggplot2 equivalents across affected functions.
- **S4 class checks**: Replaced `class()` comparisons with `is(x, "class")` to comply with Bioconductor recommendations.
- **Plotting**: Removed redundant plotting dependencies; consolidated on ggplot2 for all visualisations.
- **Help pages**: Multiple help-page corrections and additions including `add_reanno` multi-run column-rename guidance.

### Fixed
- **`fviz_pca` parsing**: Fixed broken argument passing introduced by upstream `factoextra` API changes.
- **`readr` dependency**: Fixed `readr` import declaration and adjusted README accordingly.

---

## [1.4.0] – 2025-05-28

### Added
- **`PAC_analyze` wrapper**: New high-level wrapper consolidating common analysis steps into a single call.
- **`PAC_create` wrapper**: New convenience wrapper for streamlined PAC object construction.
- **Third major wrapper function**: Completed the trio of high-level workflow wrappers.
- **`output` option in `PAC_analyze`**: Allows writing results directly to disk.

### Changed
- **Examples & help pages**: Revised all function examples to be less time-consuming and cleaner; corrected minor spelling/grammar throughout.
- **Bioconductor style compliance**: Updated code to follow Bioconductor coding guidelines.

### Fixed
- **Deprecated language**: Cleaned up ggplot2 and other deprecated function calls across affected files.
- **Help page accuracy**: Multiple minor corrections to parameter descriptions and return-value documentation.

---

## [1.1.1] – 2021-08-01

### Added
- First hard public release of `seqpac` on Bioconductor.
- S4 class compatibility throughout the package.
- `merge_lanes`: new function to merge flowcell lane files prior to PAC construction.
- `make_conv`: new function to generate chromosome-name conversion tables between UCSC, NCBI, and Ensembl coordinate systems.

### Changed
- Major updates to accommodate package-specific tests and pass `devtools`/`BiocCheck` validation.
- Streamlined PAC generation and annotation pipeline.
- Vignette updated to reflect new workflow and functions.

---

## [0.99.18] – 2021 *(Bioconductor review)*

### Changed
- Many minor updates to comply with Bioconductor reviewer feedback.

---

## [0.99.8] – 2021 *(Bioconductor review)*

### Added
- `make_count` now supports chunked, on-disk processing for low-memory/low-end systems.
- Quick-start section added to vignette.

### Changed
- More efficient function examples to reduce check times.
- Multiple minor updates addressing Bioconductor review comments.

### Fixed
- Minor bug corrections throughout.

---

## [0.99.4] – 2021 *(Bioconductor review)*

### Fixed
- Corrected notes raised by Bioconductor reviewer.
- Reduced example run times.

---

## [0.99.3] – 2021 *(Bioconductor review)*

### Fixed
- Corrected additional notes from Bioconductor reviewer.

---

## [0.99.2] – 2021 *(Bioconductor review)*

### Fixed
- Corrected errors and warnings identified by bioconductor.org automated checks.

---

## [0.99.1] – 2021 *(Bioconductor review)*

### Added
- Preparations for Bioconductor submission.
- Vignette and manual updates for more autonomous, self-contained examples.

### Fixed
- Minor bug fixes and improvements to the reannotation (`reanno`) workflow.

---

## [1.0.1] – 2020

### Added
- First GitHub release.
- Working implementation for constructing PAC objects (S3 class).
- Foundation of functions for sequence-based counting and annotation.
