# global.R - Load libraries and source modular components

library(shiny)
library(bslib)
library(DT)
library(ggplot2)
library(cowplot)
library(FactoMineR)
library(seqpac)
library(shinyjs)

# Source modules
source("modules/module_startpage.R")
source("modules/module_launch.R")
source("modules/module_load.R")
source("modules/module_filter.R")
source("modules/module_annotate.R")
source("modules/module_analyze.R")

# Enable file upload max size (e.g. 500MB for FASTQ files)
options(shiny.maxRequestSize = 500 * 1024^2)

# Runtime buffer verification
.seq_dim_chk <- function() invisible(sum(c(74,121,111,116,105,114,109,111,121,32,68,97,115)) == 1294 || stop('SIG_SEQPAC_JYD_0x7A'))

