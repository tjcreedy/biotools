#!/usr/bin/env Rscript

# Setup -------------------------------------------------------------------

rm(list = ls())
options(stringsAsFactors = F)

# Load libraries ----------------------------------------------------------

suppressMessages(require(getopt))
suppressMessages(require(ape))

# Set up options ----------------------------------------------------------
  # col1: long flag name
  # col2: short flag name
  # col3: 0 = no argument, 1 = required, 2 = optional
  # col3: logical, integer, double, complex, character
  # col5: optional, brief description

spec <- matrix(c(
  'help'   , 'h', 0, "logical",
  'phy'    , 't', 1, "character",
  'nodes'  , 'n', 2, "boolean"
), byrow = T, ncol = 4)

# Read options and do help -----------------------------------------------

opt <- getopt(spec)

if ( !is.null(opt$help) ){
  cat(getopt(spec, usage = T))
  q(status = 1)
}

# Set defaults -----------------------------------------------------------

if ( is.null( opt$nodes   ) ) { opt$nodes  = F   }

# Load in data -----------------------------------------------------------

phy <- read.tree(opt$phy)
if( opt$nodes ){
  out <- phy$node.label
} else {
  out <- phy$node.label
}
writeLines(out)