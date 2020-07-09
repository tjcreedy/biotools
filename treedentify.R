#!/usr/bin/env Rscript

# Setup -------------------------------------------------------------------

rm(list = ls())
options(stringsAsFactors = F)

# Global variables --------------------------------------------------------

taxlevels <- c("subspecies","species","superspecies",
               "subgenus","genus",
               "infratribe","subtribe","tribe","supertribe",
               "infrafamily","subfamily","family","superfamily",
               "parvorder","infraorder","suborder","order","superorder","magnorder",
               "cohort","legion",
               "parvclass","subteclass","infraclass","subclass","class","superclass",
               "microphylum","infraphylum","subphylum","phylum","superphylum",
               "infrakingdom","subkingdom","kingdom","superkingdom") # From lowest to highest

# Load libraries ----------------------------------------------------------

suppressMessages(require(getopt))
suppressMessages(require(geiger))
suppressMessages(require(taxize))
suppressMessages(require(plyr))

# Load functions ----------------------------------------------------------

listancestors <- function(tree, n, inc.n = F){
  root <- Ntip(tree)+1
  
  if(n %in% tree$tip.label){
    n <- which(tree$tip.label == n)
  }
  if(n == root){
    return(0)
  } else {
    p <- tree$edge[tree$edge[,2] == n,1]
    
    if(p == root){
      return(p)
    } else {
      out <- c(p,listancestors(tree, p, inc.n = F))
      if(inc.n == T){
        return(c(n, out))
      } else {
        return(out)
      }
    }
  }
}

# Set up options ----------------------------------------------------------
# col1: long flag name
# col2: short flag name
# col3: 0 = no argument, 1 = required, 2 = optional
# col3: logical, integer, double, complex, character
# col5: optional, brief description

spec <- matrix(c(
  'help'      , 'h', 0, "logical"   , "show this helpful message",
  'auth'      , 'a', 1, "chacter"   , "an ncbi_authentication text file with your API key as the second line not beginning with #",
  'phylo'     , 'p', 1, "character" , "the phylogeny with known and unknown taxonomy",
  'newprefix' , 'n', 2, 'character' , "all novel tips (to be identified) start with this (e.g. \'otu\')",
  'metadata'  , 'm', 2, "character" , "metadata containing phylogeny for the tree",
  'threads'   , 't', 2, 'integer'   , "number of threads to run on",
  'taxcache'  , 'c', 2, 'character' , "path to a .RDS cache of taxonomy data to read from and/or write to",
  'notstrict' , 's', 0, 'logical'   , "don't strictly taxonomise the tree"
), byrow = T, ncol = 5)

# Read options and do help -----------------------------------------------

opt <- getopt(spec)

if ( is.null(opt) | !is.null(opt$help) ){
  message(getopt(spec, usage = T))
  q(status = 1)
}

rm(spec)

# Set defaults -----------------------------------------------------------

# opt$phylo <- "3.RAxML_Mitogenome+OTUs.nwk"
# opt$auth <- "tjc_ncbi_authentication.txt"
# opt$notstrict <- TRUE
# opt$taxcach <- "ncbi_cache.RDS"

if ( is.null(opt$auth)      )  { stop("NCBI authentication file is required") }
if ( is.null(opt$phylo)     )  { stop("Phylogeny is required")                }
if ( is.null(opt$threads)   )  { opt$threads    <- 1                          }
if ( is.null(opt$notstrict) )  { opt$notstrict  <- FALSE                      }

# Load in tree and get tips ----------------------------------------------

tree <- read.tree(opt$phylo)
tips <- tree$tip

noveltips <- c()
unknowntips <- tips
if ( ! is.null(opt$newprefix) ){
  noveltips <- tips[grepl(paste0('^', opt$newprefix), tips)]
  unknowntips <- unknowntips[!unknowntips %in% noveltips]
}
rm(tips)

# Load the cache if present -----------------------------------------------

taxcache <- list()
if ( !is.null(opt$taxcache) && file.exists(opt$taxcache)) {
  taxcache <- readRDS(opt$taxcache)
}

# Load api key ------------------------------------------------------------

authlines <- readLines(opt$auth)
authlines <- authlines[! grepl("^#", authlines)]
options(ENTREZ_KEY = authlines[2])

# Parse the metadata for names and extract taxonomy -----------------------

metadataout <- c()

if ( !is.null(opt$metadata) ) { 
  # Load data
  metadata <- read.csv(opt$metadata, row.names = 1)
  
  # Extract names and match to tree
  metanames <- row.names(metadata)
  metanames <- intersect(unknowntips, metanames)
 
  # Subset only names of interest, check that we have data
  metadata <- metadata[metanames, ]
  if ( nrow(metadata) == 0 ) {
    stop("Error: cannot identify match any tip names on the tree with any values in the first column of the metadata table")
  }
  
  # Extract taxon_ids if present and store, keep remainder for taxonomy
  if ( 'taxon_id' %in% colnames(metadata) ) {
    metadataout <- metadata[, 'taxon_id']
    metadataout <- metadataout[metadataout != '' & ! is.na(metadataout)]
    metadata <- metadata[!row.names(metadata) %in% names(metadataout), ]
  } 
  
  # If any rows remain
  if ( nrow(metadata) > 0 ) {
    # Extract the present taxonomy levels
    preslevels <- taxlevels[taxlevels %in% colnames(metadata)]
    
    # Subset rows, sorting by level if present
    if ( length(preslevels) > 0 ) {
      metadata <- metadata[, preslevels]
      
      # Find taxon ids using taxize by working up the available taxonomic levels
      metadata$taxon_id <- apply ( metadata, 1, function(tax) {
        out <- NA
        for ( t in na.omit(tax[tax != '']) ) {
          out <- get_uid(t, messages = F, ask = F)[1]
          if ( !is.na(out) ) break
        }
        return(out)
      })
      
      # Drop anything that didn't return a taxon id
      metadata <- metadata[! is.na(metadata$taxon_id), 'taxon_id']
      
      # Concatenate all the found taxon_ids
      metadataout <- c(metadataout, metadata)
      
      if (nrow(metadataout) == 0) {
        stop("Error: failed to find any usable taxonomy information in the metadata file")
      }
        
    # Throw error if no taxonomy levels and we didn't previously find taxon_ids
    } else if ( length(metadataout) == 0 ) {
      stop("Error: cannot identify any columns in the metadata table matching taxon_id or known taxonomy levels")
    }
  }
  
  unknowntips <- unknowntips[! unknowntips %in% names(metadataout)]
  
  message(paste("Retrieved taxon ids for", length(metadataout), "tips from metadata"))
}

# Retrieve taxonomy based on Genbank accessions ---------------------------
gbaccessions <- c()
gbdf <- data.frame()
gbout <- c()

if ( length(unknowntips) > 0 ){
  # Extract genbank accessions from the tips
  gbaccessions <- sapply(unknowntips, function(tip){
    regmatches(tip, regexpr("(?:[^A-Za-z0-9]|^)[A-Z]{1,2}_{0,1}[0-9]{5,8}(?:[^A-Za-z0-9]|$)", tip))[1]
  })
  
  if ( sum(!is.na(gbaccessions)) > 0 ) {
    # Create dataframe and drop any rows without putative genbank accessions
    gbdf <- data.frame(gbacc = gbaccessions, id = unknowntips)[gbaccessions != '' & !is.na(gbaccessions),]
    
    # Retrieve the taxon ids from genbank and drop rows without taxon ids
    gbdf$taxon_id <- unlist(genbank2uid(gbdf$gbacc))
    gbout <- gbdf[! is.na(gbdf$taxon_id), 'taxon_id']
    names(gbout) <- row.names(gbdf)
    
    # See if we have any unknown tips left
    unknowntips <- unknowntips[! unknowntips %in% names(gbout)]
    message(paste("Retrieved taxon ids for", length(gbout), "tips from GB accession numbers in tip names"))
  }
}

rm(gbaccessions, gbdf)

# Retrieve the taxonomy directly from the tip names -----------------------
tipuids <- c()
tipout <- c()

if ( length(unknowntips) > 0 ){
  tipuids <- sapply(unknowntips, function(tip){
    out <- NA
    
    # Remove any underscores
    tiprn <- gsub('_+', ' ', tip)
    
    # Split up by common separators
    tipsplit <- strsplit(tiprn, ';,')[[1]]
    
    # If only one name still, try getting the taxon_id just by searching the whole name
    if ( length(tipsplit) == 1 ){
      out <- get_uid(tiprn, ask = F, messages = F)[1]
      
      # If no luck, try breaking it up
      if ( is.na(out) ) {
        tipsplit <- strsplit(tiprn, ' ')[[1]]
      }
    }
    
    if ( is.na(out) ) {
      # Try each subset from the end
      tipsplit <- rev(tipsplit)
      for ( t in tipsplit ) {
        out <- get_uid(t, messages = F, ask = F)[1]
        if ( !is.na(out) ) break
      }
    }
    
    return(out)
  })
  
  # Drop any that failed to return
  tipout <- tipuids[tipuids != '' & ! is.na(tipuids)]
  
  # See if we have any unknown tips left
  unknowntips <- unknowntips[! unknowntips %in% names(tipout)]
  
  message(paste("Retrieved taxon ids for", length(tipout), "tips from taxonomy in tip names"))
}

rm(tipuids)
# Bring together ----------------------------------------------------------

# Bring together unknown tips
if ( length(unknowntips) > 0 & length(noveltips) > 0 ){
  message(paste("Could not find taxon ids for", length(unknowntips), "tips, these will be added to the known novel tips"))
}
noveltips <- c(unknowntips, noveltips)
if ( length(noveltips) == 0 ){
  stop("Error: no tips to search for taxonomy for!")
}
rm(unknowntips)

# Bring together known tips
tip_taxonid <- c(metadataout, gbout, tipout)
rm(metadataout, gbout, tipout)

# Get taxonomy classification
  # Find unique
uuids <- unique(tip_taxonid)
  # Extract any from taxcache
inlocal <- uuids[uuids %in% names(taxcache)]
taxlocal <- list()
if ( length(inlocal) > 0 ) {
  taxlocal <- taxcache[inlocal]
}
  # Get from NCBI
taxncbi <- classification(uuids[! uuids %in% inlocal], db = "ncbi")
  # Concatenate
taxall <- c(taxlocal, taxncbi)
if ( ! is.null(opt$taxcache) ) {
  saveRDS(taxall, opt$taxcache)
}
  # Grab complete
taxonid_taxonomy <- taxall[tip_taxonid]
rm(taxall, taxlocal, taxncbi, uuids)

  # Convert to data frame
taxonid_taxonomy <- do.call('rbind.fill', lapply(taxonid_taxonomy, function(cls){
  data.frame(matrix(rev(cls$name[cls$rank != 'no rank']), nrow = 1, 
                    dimnames = list(cls$id[nrow(cls)], 
                                    rev(cls$rank[cls$rank != 'no rank']))))
}))
  # Sort columns by ascending level
preslevels <- taxlevels[taxlevels %in% colnames(taxonid_taxonomy)]
taxonid_taxonomy <- taxonid_taxonomy[, preslevels]
  # Bring together to full table and write out
taxonomy <- cbind(id = names(tip_taxonid), taxon_id = tip_taxonid,
                  taxonid_taxonomy)
write.csv(taxonomy, "existing_taxonomy.csv", row.names = F, quote = F)
row.names(taxonomy) <- taxonomy$id

# Apply taxonomy to phylo
message(paste("Applying taxonomy. The next message will probably be warnings from geiger::nodelabel.phylo - don't be alarmed:",
              "> The message 'redundant labels encountered' means that these levels are fully nested, no big deal.",
              "> The message 'redundant labels encountered at root' means that the entire tree falls inside these taxa, also no big deal unless the first value is a very low taxonomic level.",
              "> The message 'labels missing from phy' means that it couldn't place these taxonomic levels. We don't expect it to be able to place everything, but if this seems like a lot you could try running again with -notstrict (if you aren't already.",
              sep = '\n'))
tree_nl <- nodelabel.phylo(tree, taxonomy[preslevels], strict = !opt$notstrict, ncores = opt$threads)

write.tree(tree_nl, "taxonomised_tree.nwk")

# Get taxonomy for each novel tip
noveltaxonomy <- lapply(noveltips, function(tip){
  # Get ancestor nodes for this tip
  ancnods <- listancestors(tree_nl, tip)
  ancnods <- rev(sort(ancnods))
  # Retrieve any names
  nodes <- tree_nl$node.label[ancnods - length(tree$tip.label)]
  nodes <- nodes[nodes != '']
  out <- NA
  if ( length(nodes) > 0 ) {
    est <- F
    # Take the lowest taxonomic level
    lowest <- nodes[1]
    # Check if a guessed node
    if ( grepl('\"', lowest) ) {
      est <- T
      lowest <- gsub('\"', '', lowest)
    }
    taxid <- get_uid(lowest, ask = F, messages = F)[1]
    # Find the taxonomic level this name belongs to
    j <- which(apply(taxonomy, 2, function(names) lowest %in% names))[1]
    # Return the fullest available taxonomy
    out <- cbind(taxon_id = taxid, est = est, taxonomy[taxonomy[,j] == lowest & !is.na(taxonomy[,j]), j:ncol(taxonomy)][1, ])
    row.names(out) <- NULL
    
  }
  return(out)
})

success <-  sapply(noveltaxonomy, is.data.frame)
noveltaxonomy <- do.call('rbind.fill', noveltaxonomy[success])
row.names(noveltaxonomy) <- noveltips[success]
preslevels <- taxlevels[taxlevels %in% colnames(noveltaxonomy)]
noveltaxonomy <- noveltaxonomy[, c('taxon_id', 'est', preslevels)]

message(paste("Successfully assigned taxonomy to", sum(success), "out of", length(success), "novel or unknown tips."))
write.csv(noveltaxonomy, "new_taxonomy.csv", row.names = T, quote = F)

# Plot tree ---------------------------------------------------------------

# alltax <- rbind.fill(taxonomy, noveltaxonomy)
# row.names(alltax) <- c(row.names(taxonomy), row.names(noveltaxonomy))
# alltax[is.na(alltax)] <- ''
# alltax$label <- paste0(row.names(alltax), '->', apply(alltax[,c('class', 'order', 'family', 'genus')], 1, function(x) paste(x[x!=''], collapse = ';')), 
#                       ifelse(alltax$est == 'TRUE', '?', ''))
# alltax$label <- gsub('Insecta;', '', alltax$label)
# 
# tree_nl$tip.label <- alltax[tree_nl$tip.label, 'label']
# plot.phylo(ladderize(tree_nl,right=FALSE), cex=0.35, 
#            type = 'fan', label.offset = 2.5,
#            no.margin=TRUE, edge.color="gray", edge.width=0.5)
# nodelabels(pie = ifelse(tree_nl$node.label != '', 1, NA), cex = 0.1)
# nodelabels(tree_nl$node.label, cex=0.45, adj = c(1.2, 0.5), col="red", frame="n")

q('no', 0, F)