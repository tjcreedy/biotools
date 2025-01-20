# Setup -------------------------------------------------------------------

rm(list = ls())
options(stringsAsFactors = F)

# Load libraries ----------------------------------------------------------

suppressMessages(require(getopt))
suppressMessages(require(tidyr))
suppressMessages(require(dplyr))
suppressMessages(require(magrittr))
suppressMessages(require(purrr))
suppressMessages(require(ape))
suppressMessages(require(tibble))
suppressMessages(require(knitr))

# Load functions ----------------------------------------------------------

# From phangorn
upgma <- function (D, method = "average", ...) 
{
  DD <- as.dist(D)
  hc <- hclust(DD, method = method, ...)
  result <- as.phylo(hc)
  result <- reorder(result, "postorder")
  result
}

cluster <- function(.x, .y, distm){
  taxa <- .y %>% unlist() %>% {.[!is.na(.)]} %>% rev()
  taxon <- taxa  %>% pluck(1)
  if(nrow(.x) > 1){
    d <- dist.dna(ASVs[.x$ASV_ID], distm)
    disttaxon <- NA
    if( is.null(taxon) ){
      dists <- ddist
      disttaxon <- c("[no taxonomy]" = "[no taxonomy]")
      taxon <- "notax"
    } else {
      for(i in 1:length(taxa)){
        #i = 1
        if(any(similaritylin[, names(taxa)[i]] == taxa[i] & !is.na(similaritylin[, names(taxa)[i]]))){
          dists <- left_join(.y, similaritylin, by = c(names(taxa[i:length(taxa)]))) %>% pull(maxdist)
          disttaxon <- taxa[i]
          break
        }
      }
    }
    cutree(as.hclust.phylo(upgma(d)), h = mean(dists)) %>% enframe() %>%
      setNames(c("ASV_ID", "otu")) %>%
      mutate(otu = paste0(taxon, " otu", sprintf("%03d", otu)),
             distance = mean(dists), disttaxon = disttaxon, distrank = names(disttaxon)) %>%
      return
  } else {
    .x %>% ungroup %>% select(ASV_ID) %>% mutate(otu = paste0(taxon, " otu001"), 
                                                 distance = NA, disttaxon = NA, distrank = NA) %>% return
  }
}

# Set global variables ----------------------------------------------------

taxlevels <- readLines("https://raw.githubusercontent.com/tjcreedy/constants/main/taxlevels.txt")

# Set up options ----------------------------------------------------------
# col1: long flag name
# col2: short flag name
# col3: 0 = no argument, 1 = required, 2 = optional
# col3: logical, integer, double, complex, character
# col5: optional, brief description

spec <- matrix(c(
  'help'        , 'h', 0, "logical"  , "show this helpful message.",
  'asvfasta'    , 'a', 1, "character", "path to a fasta containing the sequences for each asv.",
  'taxonomy'    , 't', 1, "character", "path to a csv table of taxonomy for each asv, the first column (ASV_ID or qseqid) must contain asv ids and all other columns taxonomic levels (species, genus, etc). Missing data should be blank", 
  'reads'       , 'r', 1, "character", "path to a csv table recording the read count for each asv within each sample, the first column (ASV_ID) must contain asv names and all subsequent columns samples.",
  'distmodel'   , 'm', 2, "character", "the evolutionary model for sequence dissimilarity computation, see the model options here https://www.rdocumentation.org/packages/ape/versions/5.8/topics/dist.dna, default = \"raw\"",
  'defaultdist' , 'd', 2, "double"   , "the distance at which to cluster ASVs without taxonomy. If not supplied, ASVs without taxonomy will be removed (recommended)",
  'outprefix'   , 'o', 2, "character", "prefix file path to write clustering data."
), byrow = T, ncol = 5)

# Read options and do help -----------------------------------------------

opt <- getopt(spec)

if ( is.null(opt) | !is.null(opt$help) ){
  message(getopt(spec, usage = T))
  q(status = 1)
}

rm(spec)

# opt$taxonomy <- "filtered_taxonomy.csv"
# opt$reads <- "filtered_readstable.csv"
# opt$asvfasta <- "ampliseqout_2024-02-20/codon_filter/ASV_codon_filtered.fna"
# opt$distmodel <- "raw"
# opt$outprefix <- "./"
# opt$defaultdist <- 0.03

if( is.null(opt$distmodel) ) opt$distmodel <- "raw"
if( is.null(opt$outprefix) ) opt$outprefix <- "./"

# Load data and check -----------------------------------------------------

taxonomy <- read.csv(opt$taxonomy)
reads <- read.csv(opt$reads)
ASVs <- read.FASTA(opt$asvfasta)

# Check colnames

if(names(taxonomy)[1] == "qseqid"){
  names(taxonomy)[1] <- "ASV_ID"
} else if(names(taxonomy)[1] != "ASV_ID"){
  stop(" the first column of the supplied taxonomy csv is not called \"ASV_ID\" or \"qseqid\".")
}

if(names(reads)[1] != "ASV_ID"){
  stop(" the first column of the supplied reads csv is not called \"ASV_ID\".")
}

present_taxlevels <- names(taxonomy)[-1]
taxnamerr <- present_taxlevels[!present_taxlevels %in% taxlevels]
if(length(taxnamerr) > 0){
  message(paste("\nWarning: some column names in the supplied taxonomy csv are not known taxonomic levels. Taxonomic levels must be lower case. The following columns will be dropped from the outputs. Issues:",
             paste(taxnamerr, collapse = ",")))
}

if(! "species" %in% present_taxlevels ){
  stop(" species is not present in the supplied taxonomy")
}

taxlevels <- taxlevels[taxlevels %in% present_taxlevels]

# Check that all ASV IDs match up
if(! ( all(taxonomy$ASV_ID %in% names(ASVs)) & all(taxonomy$ASV_ID %in% reads$ASV_ID) )  ){
  stop(" one or more ASVs in the supplied taxonomy csv are missing from the reads file and/or ASV fasta. Please check")
}


# Convert blanks to NAs
taxonomy %<>% 
  select(c(ASV_ID, all_of(taxlevels))) %>%
  mutate(across(!ASV_ID, ~ na_if(.x, "")))

# Create long version
taxonomylong <- taxonomy %>% 
  pivot_longer(!ASV_ID, names_to = "rank", values_to = "clade") %>%
  mutate(rank = factor(rank, levels = taxlevels)) %>%
  filter(!is.na(clade))

message("\nSuccessfully read in taxonomy data for ", nrow(taxonomy), " ASVs and retrieved sequences.")

# Find unique taxonomies ----------------------------------------------------------------------

lineages <- taxonomy %>%
  select(-ASV_ID) %>%
  unique

# How many assigned to different taxonomic levels?
message("\nInput taxonomy has the following numbers of ASVs assigned at each rank:")
rankcounts <- taxonomylong %>% count(rank) 
rankcounts %>% kable %>% print

if( ! "species" %in% rankcounts$rank ){
  stop(" no ASVs are assigned to species level, cannot perform clustering.")
}

notax <- taxonomy$ASV_ID[! taxonomy$ASV_ID %in% taxonomylong$ASV_ID]
if( length(notax) > 0 ){
  msg <- paste("\nWarning:", length(notax), "ASVs have no taxonomy, these will be")
  if( is.null(opt$defaultdist) ){
    msg <- paste(msg, "removed (recommended). To retain and cluster these ASVs, set -d/--defaultdist.")
    taxonomy <- taxonomy %>% 
      filter(!ASV_ID %in% notax)
  } else {
    msg <- paste(msg, "clustered using a distance of", paste0(opt$defaultdist, "."))
  }
  message(msg)
}

reads %<>% filter(ASV_ID %in% taxonomy$ASV_ID)
ASVs <- ASVs[names(ASVs) %in% taxonomy$ASV_ID]

if( rankcounts %>% filter(rank == "species") %>% pull(n) == nrow(taxonomy) ){
  message("\nAll ASVs are assigned to species, so skipping similarity clustering and just grouping ASVs by species.")
  taxonomy_out <- taxonomy %>%
    mutate(otu_final = species)
  
} else {

# Find similarity for ASVs assigned to the same species ---------------------------------------
message("\nCalculating within-species similarity.")
similarity <- taxonomylong %>%
  filter(rank == "species") %>%
  group_by(clade) %>%
  mutate(n = n()) %>%
  filter(n > 1) %>%
  group_map( ~ data.frame(clade = .y, maxdist = max(dist.dna(ASVs[.x$ASV_ID], opt$distmodel)))) %>%
  list_rbind()

similaritylin <- 
  right_join(lineages, similarity, by = c("species" = "clade")) %>%
  arrange(maxdist)

# How many species to compute maximum within-species distance
message(nrow(similarity), " species have >1 ASV and were used for calculating maximum within-species similarity.")

# Find ASVs not assigned to species and cluster them ------------------------------------------
message("\nPerforming clustering...")
groups <- taxlevels[taxlevels != "species"] %>% rev

clusters <- taxonomy %>%
  filter(is.na(species)) %>%
  select(-species) %>%
  group_by(pick(all_of(groups))) %>% 
  group_map(cluster, distm = opt$distmodel, ddist = opt$defaultdist) 

# DEBUG
# tocluster <- taxonomy %>% filter(is.na(species)) %>%  select(-species) %>% group_by(pick(all_of(groups))) %>% group_map( ~list(".x" = .x, ".y" = .y))
# for(i in 1:length(tocluster)){
#     cluster(tocluster[[i]]$.x, tocluster[[i]]$.y, opt$distmodel, opt$defaultdist)
# }
# i

distdata <- map(clusters, ~slice(.x, 1) %>% select(starts_with("dist"))) %>% 
  bind_rows %>% 
  filter(!is.na(distance)) %>%
  mutate(distrank = factor(distrank, levels = taxlevels))

clusters %<>%
  list_rbind() %>%
  rename_with( ~ paste0("otu_", .x), starts_with("dist"))

# Check
# clusters %>% filter(otu_disttaxon == "[no taxonomy]")

# Analyse clustering ------------------------------------------------------

# How many ASVs clustered?
message("\n\nClustered ", nrow(clusters), " ASVs not assigned to species into ", 
        clusters %>% pull(otu) %>% unique %>% length, " OTUs")

# How many ASVs clustered automatically because single within taxon
nsingletaxa <- clusters %>% filter(is.na(otu_distance)) %>% count(otu) %>% nrow
message("Of these, ", nsingletaxa, " ASVs were the only representative of their lowest taxon, so formed single-ASV OTUs")

# How many groups of clustering performed
message("The remaining ", clusters %>% filter(!is.na(otu_distance)) %>% nrow(), 
        " ASVs were clustered into ", nrow(distdata), " OTUs")

# How many groups used which rank?
message("In forming these ", nrow(distdata), " OTUs, clustering was performed based on closest relatives at the following levels:")
distdata %>% count(distrank) %>% kable %>% print

message("\nSummary of within-species maximum distances used for clustering:")
print(summary(distdata$distance))

# Assign clustering to input data -----------------------------------------

taxonomy_out <- taxonomy %>%
  left_join(clusters, by = "ASV_ID") %>%
  mutate(otu_final = ifelse(is.na(otu), species, otu))

}

# Generate outputs --------------------------------------------------------

taxonomy_out %>%
  select(otu_final, all_of(taxlevels)) %>%
  unique() %>% 
  write.csv(paste0(opt$outprefix, "taxonomy_collapsed.csv"), row.names = F)
write.csv(taxonomy_out, paste0(opt$outprefix, "taxonomy_full.csv"), row.names = F)

reads_out <- reads %>%
  left_join(taxonomy_out %>% select(ASV_ID, otu_final), by = "ASV_ID")

reads_out %>%
  select(-ASV_ID) %>%
  relocate(otu_final) %>%
  group_by(otu_final) %>%
  summarise(across(everything(), sum), .groups = "drop") %>%
  write.csv(paste0(opt$outprefix, "reads_collapsed.csv"), row.names = F)
write.csv(reads_out %>% relocate(ASV_ID, otu_final), paste0(opt$outprefix, "reads_full.csv"), row.names = F)

seqrename <- reads_out %>%
  pivot_longer(!c(ASV_ID, otu_final)) %>%
  group_by(otu_final, ASV_ID) %>%
  summarise(reads = sum(value), .groups = "drop") %>%
  group_by(otu_final) %>%
  filter(reads == max(reads)) %>%
  {setNames(pull(., otu_final), pull(., ASV_ID))}

outASVs <- ASVs[names(seqrename)]
names(outASVs) <- unname(seqrename[names(outASVs)])

write.FASTA(outASVs, paste0(opt$outprefix, "otu_sequences.fasta"))

message("\nInput reads and taxonomies output as-is with OTU assignments appended to _full csv files. Taxonomy and reads collapsed to OTUs output to _collapsed csv files. OTU sequences output. All done!")
