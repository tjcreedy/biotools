
# Introduction ------------------------------------------------------------

# By Thomas J. Creedy, 2020-04-03

# This script is a template for analysing metabarcoding data generated
# using the current Vogler Lab metabarcoding pipeline.

# I strongly recommend using these commands as a guide only and ensuring
# you understand what each command does and how to modify it to answer your
# specific questions.

# Not everything here will work for all data. Please try and figure out how
# to adapt these commands to your data - the purpose of this script is not
# to be a comprehensive universal analysis tool but to be a template to get
# you started.

# Setup -------------------------------------------------------------------
  # Commands in this section clean up from any previous runs and implements
  # any global settings.

  # Remove all objects from the environment
rm(list = ls())

  # Set global options
options(stringsAsFactors = F)


# Load libraries ----------------------------------------------------------
  # Load the libraries (packages) you need for your data manipulation and 
  # analysis. Remember, the order of loading is sometimes important.

library(plyr)
library(dplyr)
library(reshape2)
library(vegan)
library(ggplot2)

# Load functions ----------------------------------------------------------
  # Load any custom functions. These are assumed to be in the root of your
  # working directory.

source("check_expected_richness.R")
source("rarexplore.R")

# Load data ---------------------------------------------------------------
  # Load the three core data tables. If you have more than these data,
  # great! Add loading commands for them here. I would generally suggest
  # merging and/or reconfiguring them as necessary to be similar to the
  # configuration of these tables.


  # Load the reads table output by vsearch --usearch_global --tabbedout

reads <- read.table("reads.tsv", header = T, sep = "\t", row.names = 1, comment.char = '')

  # Load metadata table (change command if yours is a different format!)
    # Note that here I'm assuming that the first column of the metadata
    # table is the sample names.

metadata <- read.csv("metadata.csv", row.names = 1)

    # If the sample names are in a different column, do this:

#metadata <- read.csv("metadata.csv")
#row.names(metadata) <- metadata$samplenames
#metadata <- subset(metadata, select = -samplenames)

  # Load taxonomy table
    
  # From SINTAX
taxonomy <- read.table("taxonomy.tsv", sep = "\t", row.names = 1)
colnames(taxonomy) <- c("taxonomy", "strand", "selected")


# Organise data -----------------------------------------------------------
  # Here we make sure all the data corresponds properly to each other and 
  # any reconfiguration or merging is done as needed.

  # Generally, samples should be rows and species should be columns, so
  # we transpose the usearch_global --tabbedout reads table

reads <- t(reads)

  # If we have multiple metabarcoding samples that actually represent the 
  # same ecological sample, here we might merge them together

## Some merging using rowsum(reads, group = XXX)

  # Next, check that our metadata and reads table correspond. For both 
  # tables, the sample names are in the row.names
    
    # Are all MBC samples in metadata? If not, problem!!!
row.names(reads) %in% row.names(metadata) 

    # Make the metadata correspond to the reads
metadata <- metadata[row.names(reads), ]
      # This keeps only metadata rows where the sample occurs in MBC, and 
      # sorts the rows to match the MBC row order

  # The data in the taxonomy table needs to be filtered and separated out
    # As above
colnames(reads) %in% row.names(taxonomy)
taxonomy <- taxonomy[colnames(reads),]

    # We create a detailed version of the taxonomy table with the scores
    # This code is based on Yige Sun's rewriting of this step, thanks!
      # Add the OTUs as a column
taxonomy$otu <- row.names(taxonomy)
      # Separate the taxonomy column into multiple columns by the ',' character
taxdetailed <- separate(taxonomy, taxonomy, paste0('V', 1:(max(str_count(taxonomy$taxonomy, ','))+1)), 
                        sep = ",", fill = 'right', remove = T) %>% 
      # Turn the separate columns for taxa into rows
  pivot_longer(., cols = starts_with('V'), values_to = 'taxon', values_drop_na = T) %>%
      # Retain only the OTU and taxon columns
  select(otu, taxon) %>%
      # Separate the taxon column into three columns by colon and parentheses
  separate(., taxon, into = c('level', 'taxon', 'score', NA), sep = "[:()]")

    # Record which taxa were selected
taxdetailed$selected <- taxdetailed$score == 1

    # Rename the taxonomic levels
levelnames <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
names(levelnames) <- substr(levelnames, 1, 1)
taxdetailed$level <- factor(levelnames[taxdetailed$level], levels = levelnames)

    # Now we overwrite the original taxonomy table with the detailed one
taxonomy <- dcast(taxdetailed, otu ~ level, value.var = "taxon")

    # We may only wish to include taxonomy info if it's above a certain score
#taxonomy <- dcast(taxdetailed[taxdetailed$score > 0.65, ], otu ~ level, value.var = "taxon", fill = '')
    
    # Set the row names
row.names(taxonomy) <- taxonomy$otu
taxonomy <- subset(taxonomy, select = -otu)    

    # Finally, re-sort the taxonomy
taxonomy <- taxonomy[colnames(reads), ]

rm(taxdetailed)

# Filtering data ----------------------------------------------------------
  # This is the most important step, where you remove irrelevant or likely
  # incorrect data. Depending on your questions you could apply none, some
  # or all of these steps.

  # Taxonomic filtering - retain only OTUs of a certain taxon, e.g.:

taxonomy <- taxonomy[taxonomy$phylum == "Arthropoda", ]
taxonomy <- taxonomy[taxonomy$order %in% c("Coleoptera", "Hemiptera"), ]

reads <- reads[, rownames(taxonomy)]

  # Read number filtering 
    # Set counts less than a certain value to 0 - this isn't recommended
    # because it doesn't take into account total reads per sample. It also
    # generates difficulties with some of the accumulation tests later

reads[reads < 2] <- 0
    
  # Read proportion filtering
    # Set counts that are less than a certain proportion of reads in a 
    # sample to 0. The threshold here is 0.5%

threshold <- 0.005

reads <- apply(reads, 1, function(s){
  s[s/sum(s) < threshold] <- 0
  return(s)
}) %>% t()
rm(threshold)


  # After filtering you may have some OTUs or samples that now have 0 reads
  # so these should be filtered out and the relevant lines in the metadata/
  # taxonomy dropped as well

reads <- reads[rowSums(reads) > 0, colSums(reads) > 0]
metadata <- metadata[rownames(reads), ]
taxonomy <- taxonomy[colnames(reads), ]


#  Standardisation --------------------------------------------------------
  # It is crucial with metabarcoding that you ensure your samples are 
  # comparable, otherwise your analyses may be invalid. There are two ways
  # we can consider whether samples are standard:
    # They have received equivalent sampling depth
    # They have accumulated a representative sample of a community
  # These assumptions may vary depending on your question

  # Check the read numbers - have all samples been equally sequenced?
    # Probably not!
hist(rowSums(reads))

  # Check the accumulation - have all samples recovered a good proportion
  # of their expected OTU richness?
    # Here we set the threshold proportion of recovery to 85%. This is 
    # probably the lowest this should be.
check <- check_expected_richness(reads, 0.85)

  # For most community ecology, the best route is to remove samples that 
  # did not reach the threshold. check_expected_richness returns a vector
  # of samples that passed this as the first item
reads <- reads[check[[1]], ]

  # An alternative to this standardisation is to rarefy. This process
  # removes reads from samples until it reaches a threshold. This 
  # standardises sampling depth but at the loss of real community data.
  # It is more suitable for methodological analyses and is not recommended
  # for community ecology.
  
    # First lets review how many samples/OTUs might be lost at different 
    # rarefaction target values:
target_values = round(seq(min(rowSums(reads)), max(rowSums(reads)), length.out = 30))
rarexplore(reads, target_values)

    # Pick a value that is a good trade-off of remaining samples and
    # remaining OTUs, say 2500, and drop samples with below this value
reads <- reads[rowSums(reads) >= 2500, ]
reads <- rrarefy(reads, 2500)

  # The most crucial part of standardisation is to convert your data into
  # presence/absence. Read numbers are not true abundances and are rarely
  # comparable between OTUs and are only comparable between samples if
  # rarefaction has been done (which is rarely appropriate). The only time
  # when you may retain read numbers is if you are interested in looking 
  # at variation between one or a small number of similar OTUs in rarefied 
  # samples.
reads[reads > 0] <- 1

  # Finally, we need to clean up after standardisation
    # Remove any OTUs which no longer have any reads
reads <- reads[, colSums(reads) > 0]

    # Remove any dropped samples/OTUs from metadata/taxonomy tables
metadata <- metadata[rownames(reads), ]
taxonomy <- taxonomy[colnames(reads), ]

