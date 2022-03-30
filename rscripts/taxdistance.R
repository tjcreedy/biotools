#!/usr/bin/env Rscript

# Setup -------------------------------------------------------------------

rm(list = ls())
options(stringsAsFactors = F)

# Load libraries ----------------------------------------------------------

suppressMessages(require(getopt))
suppressMessages(require(taxize))
suppressMessages(require(plyr))
suppressMessages(require(zoo))

# Global variables ----------------------------------------------------------------------------

# Get levels
taxlevels <- readLines("https://raw.githubusercontent.com/tjcreedy/constants/main/taxlevels.txt")
mainlevels <- c("species", "genus", "family", "order", "class", "phylum", "kingdom", "domain")

# Set up distances
maindists <- setNames(1:length(mainlevels), mainlevels)
taxdists <- setNames(na.approx(c(0, maindists[taxlevels]))[-1], taxlevels)
taxdists <- taxdists/(max(maindists)+1)

rm(mainlevels, maindists)

# Load functions ----------------------------------------------------------

validate_taxcache <- function(taxcache, path = NULL, exit.text = NULL){
  baseuids <- sapply(taxcache, function(tc) rev(tc$id)[1])
  check <- names(taxcache) == baseuids
  if( ! all(check) ){
    taxcache <- taxcache[check]
    if( ! is.null(path) ){
      saveRDS(taxcache, path)
    }
    if( ! is.null(exit.text) ){
      stop(exit.text)
    }
  }
  return(taxcache)
}

update_taxcache <- function(taxids, taxcache, auth, indent = ""){
  taxids <- taxids[! taxids %in% names(taxcache) ]
  if( length(taxids) == 0 ){
    return(list())
  }
  message(indent, "Updating taxonomy cache for new ", length(taxids), " taxids")
  newtaxcache <- get_taxonomy_from_taxids(taxids, taxcache, auth, indent = paste0(indent, "\t"))[[2]][taxids]
  return(newtaxcache)
}


get_taxonomy_from_taxids <- function(taxids, taxcache, auth, indent = ""){
  taxids <- as.character(taxids)
  uuids <- unique(taxids)
  
  # Extract any from taxcache
  taxlocal <- list()
  inlocal <- uuids[uuids %in% names(taxcache)]
  if ( length(inlocal) > 0){
    taxlocal <- taxcache[inlocal]
    message(indent, "Found taxonomy from local NCBI cache for ", length(taxlocal), " unique taxids")
  }
  
  # Get from NCBI
  uuids <- uuids[! uuids %in% names(taxlocal)]
  taxncbi <- list()
  if ( length(uuids) > 0 ){
    message(indent, "Running remote NCBI search on ", length(uuids), " taxids...", appendLF = F)
    start <- Sys.time()
    suppressWarnings(suppressMessages(
      taxncbi <- classification(uuids, db = "ncbi")
    ))
    stop <- Sys.time()
    taxncbi <- taxncbi[!is.na(taxncbi)]
    searchtime <- as.numeric(difftime(stop, start, units = "secs"))
    message("completed in ", round(searchtime, 1), " seconds, found taxonomy for ", length(taxncbi), " taxids")
    if ( is.null(auth) ) Sys.sleep(0.5)
  }
  
  # Concatenate
  taxall <- c(taxlocal, taxncbi)
  taxall <- setNames(taxall[as.character(taxids)], names(taxids))
  
  # Convert to data frame
  taxonomy <- do.call('rbind.fill', lapply(taxall, function(cls){
    data.frame(matrix(rev(cls$name[cls$rank != 'no rank']), nrow = 1, 
                      dimnames = list(cls$id[nrow(cls)], 
                                      rev(cls$rank[cls$rank != 'no rank']))))
  }))
  
  # Order data frame
  taxonomy <- taxonomy[, taxlevels[taxlevels %in% colnames(taxonomy)]]
  
  return(list(taxonomy, c(taxcache, taxncbi)))
}


# Set up options ----------------------------------------------------------
# col1: long flag name
# col2: short flag name
# col3: 0 = no argument, 1 = required, 2 = optional
# col3: logical, integer, double, complex, character
# col5: optional, brief description

spec <- matrix(c(
  'help'        , 'h', 0, "logical"  , "show this helpful message",
  'taxonomy'    , 'y', 1, "character", "path to tab-delimited table, either a) 2-columns without 
                                        column headers comprising id and NCBI taxid in that order 
                                        or b) with column headers comprising id and taxonomic ranks",
  'outputtable' , 't', 2, "character", "path to write output distance matrix as table",
  'outputdist'  , 'd', 2, "character", "path to write output distance matrix as dist object in RDS",
  'taxout'      , 'x', 2, "character", "optional, path to write out full taxonomy if input taxonomy
                                        is NCBI taxids",
  'taxcache'    , 'c', 2, 'character', "path to a RDS cache of taxonomy data to read from and/or 
                                        write to",
  'auth'        , 'a', 2, "character", "an ncbi_authentication text file with your API key as the 
                                        second line not beginning with #"
), byrow = T, ncol = 5)
spec[,5] <- gsub("\n[\t ]*", "", spec[,5])

# Read options and do help -----------------------------------------------

opt <- getopt(spec)

if ( is.null(opt) | !is.null(opt$help) ){
  message(getopt(spec, usage = T))
  q(status = 1)
}

rm(spec)

# Read and check taxonomy ---------------------------------------------------------------------

tax <- read.table(opt$taxonomy)
usetaxid <- ncol(tax) == 2
if(usetaxid){
  if(! is.numeric(tax[,2])){
    stop("Error: read two-column taxonomy table, assuming id and NCBI taxid, but second column is not numeric!")
  }
   message("Read two-column taxonomy table, assuming id and NCBI taxid")
   names(tax) <- c("id", "taxid")
} else {
  if( ncol(tax) == 1){
    stop("Error: taxonomy table only has one column, or is not tab-delimited")
  }
  message("Read multi-column taxonomy table, assuming id and taxonomic levels")
  colistaxlevel <- names(tax) %in% taxlevels
  if( sum(!colistaxlevel) != 1){
    stop("Error: more than one column in taxonomy table does not match a known taxonomic level: ",
         paste(names(tax)[!colistaxlevel]), collapse = ",")
  }
  tax <- cbind(tax[,!colistaxlevel], tax[,na.omit(match(taxlevels, names(tax)))])
  names(tax)[1] <- "id"
  tax[tax == ""] <- NA
}

# Set defaults and error check -----------------------------------------

if(usetaxid){
  if(is.null(opt$auth)){
    stop("Error: NCBI authentication is required if supplying taxids")
  }
} 

if( is.null(opt$outputtable) & is.null(opt$outputdist) ){
  stop("Error: specify at least one of --outputtable and/or --outputdist")
}

# Get taxonomy from taxid if needed -----------------------------------------------------------

if(usetaxid){
  
  # Get cache if present 
  taxcache <- list()
  if ( !is.null(opt$taxcache) && file.exists(opt$taxcache) ){
    taxcache <- readRDS(opt$taxcache)
    taxcache <- validate_taxcache(taxcache, opt$taxcache)
    message("Read supplied taxonomy cache")
  }
  
  # Load API key
  authlines <- readLines(opt$auth)
  authlines <- authlines[! grepl("^#", authlines)]
  options(ENTREZ_KEY = authlines[2])
  rm(authlines)
  
  # Get taxonomy
  gtftreturn <- get_taxonomy_from_taxids(tax$taxid, taxcache, opt$auth)
  tax <- data.frame(id = tax$id, gtftreturn[[1]])
  
  # Update cache if present
  if ( ! is.null(opt$taxcache) ) {
    newcache <- gtftreturn[[2]]
    newcache <- newcache[!names(newcache) %in% names(taxcache)]
    taxcache <- c(taxcache, newcache)
    saveRDS(taxcache, opt$taxcache)
  }
  
  rm(taxcache, gtftreturn, newcache)
  
  # Write out taxonomy if requested
  if( !is.null(opt$taxout) ){
    write.table(tax, opt$taxout)
  }
}

# Calculate distances -------------------------------------------------------------------------

makedist <- function(x, size = NULL, labels = NULL, method = NULL){
  if( is.null(size) & is.null(labels) ){
    stop("requires at least one of size or labels")
  } else if ( ! is.null(labels) ){
    attr(x, "Size") <- length(labels)
    attr(x, "Labels") <- labels
  } else if ( ! is.null(size) ){
    attr(x, "Size") <- size
  }
  attr(x, "Diag") <- attr(x, "Upper") <- F
  if( ! is.null(method) ){
    attr(x, "method") <- method
  }
  class(x) <- "dist"
  return(x)
}

taxdist <- function(taxonomy, taxdists){
  pt <- names(taxonomy)
  taxonomy <- as.matrix(taxonomy)
  # Work through all taxonomies
  rawdist <- setNames(lapply(1:(nrow(taxonomy)-1), function(i){
    taxi <- taxonomy[i, ]
    # Subset taxonomy to only those cells relevant for comparison to this record
    taxcompare <- taxonomy[i:nrow(taxonomy), !is.na(taxi), drop = F]
    taxi <- taxcompare[1,]
    out <- NULL
    if(ncol(taxcompare) == 0){
      # If no columns, distance is maximumdistance
      out <- rep(1, nrow(taxcompare)-1)
    } else if(ncol(taxcompare) == 1){
      # If only one column, can do a shortcut for calculating distance
      # If doesn't match here, the non-matches are no hit so should have maximum distance
      out <- ifelse(taxcompare[-1,1] == na.omit(taxi), taxdists[colnames(taxcompare)], 1)
    } else {
      out <- apply(taxcompare[-1,,drop = F], 1, function(taxj){
        matches <- taxi == taxj & !is.na(taxj)
        return(taxdists[colnames(taxcompare)[which(matches)[1]]])
      })
    }
    return(out)
  }), rownames(taxonomy)[-nrow(taxonomy)])
  return(makedist(unlist(rawdist), labels = rownames(taxonomy)))
}

rownames(tax) <- tax$id
distances <- taxdist(tax[, -1], taxdists)


# Output --------------------------------------------------------------------------------------

opt$outputdist = "7_dupnt_all_blastn-ex-distance.RDS"
if( !is.null(opt$outputdist) ){
  saveRDS(distances, opt$outputdist)
}
