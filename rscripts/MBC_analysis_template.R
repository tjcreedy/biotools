
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

library(tidyr)
library(plyr)
library(dplyr)
library(reshape2)
library(stringr)
library(vegan)
library(ggplot2)
library(ape)

# Load functions ----------------------------------------------------------
  # Load custom functions.

rarexplore <- function(map, size, seed = 12345){
  require(ggplot2)
  require(reshape2)
  require(plyr)
  require(magrittr)
  require(scales)
  require(vegan)
  raresample = size[2]
  rared <- lapply(size, function(raresample){
    set.seed(seed)
    raresample <- round(raresample)
    rared <- suppressWarnings(rrarefy(map, raresample))
    passrows <- rowSums(rared) >= raresample
    passcols <- colSums(rared) > 0
    rared <- rared[passrows, passcols] %>%
      matrix(., nrow = sum(passrows), ncol = sum(passcols),
             dimnames = list(names(passrows)[passrows], names(passcols)[passcols]))
    return(setNames(c(raresample,dim(rared)),
                    c("raresamples","Samples","OTUs"))) 
  }) %>% do.call(rbind, .) %>%
    data.frame() %>%
    melt(., id.vars = "raresamples", variable.name = "measure", value.name = "value")
  
  rared$value <- ifelse(rared$measure == "Samples", rared$value/nrow(map), rared$value/ncol(map))
  
  ggplot(data = rared, aes(x = raresamples, y = value))+
    geom_point()+
    theme_bw()+ 
    geom_line()+ 
    scale_y_continuous(labels = scales::percent)+
    labs(x = "Rarefaction sample", y = "Percentage of total remaining")+
    facet_wrap(~measure, scales = "free_y")
  
}

check_expected_richness <- function(map, threshold, method = "chao1", bins = 30){
  require(vegan)
  require(ggplot2)
  require(ggthemes)
  require(cowplot)
  
  stats <- data.frame(cbind(rowSums(map), t(estimateR(map))))
  colnames(stats)[1] <- "reads"
  stats$p_complete <- stats$S.obs/stats[,paste0("S.",method)]
  stats$p_complete_min <- stats$S.obs/(stats[,paste0("S.",method)] + stats[,paste0("se.", method)])
  stats$p_complete_max <- stats$S.obs/(stats[,paste0("S.",method)] - stats[,paste0("se.", method)])
  stats$p_complete_max[stats$p_complete_max > 1] <- 1
  
  passing <- data.frame(x = c(threshold - 0.05, threshold + 0.05),
                        y = c(0.5, 0.5),
                        label = c(paste0(round(sum(stats$p_complete < threshold, na.rm = T)/nrow(stats)*100), "%"),
                                  paste0(round(sum(stats$p_complete > threshold, na.rm = T)/nrow(stats)*100), "%")),
                        col = LETTERS[1:2])
  
  plots1 <- list(ggplot(stats, aes(x = p_complete))+
                   labs(y = "Count of samples") +
                   geom_histogram(bins = bins),
                 ggplot(stats)+
                   labs(y = "Cumulative proportion of samples <= x")+
                   stat_ecdf(aes(x = p_complete))+
                   stat_ecdf(aes(x = p_complete_min), linetype = 2) + 
                   stat_ecdf(aes(x = p_complete_max), linetype = 2) + 
                   geom_text(data = passing, aes(x = x, y = y, label = label, col = col))+
                   scale_colour_manual(values = c("red", "green"), guide = F))
  plots2 <- list(ggplot(stats, aes(x = reads, y = p_complete, ymin = p_complete_min, ymax = p_complete_max))+
                   labs(x = "Number of reads") +
                   scale_x_continuous(trans = 'sqrt') +
                   geom_errorbar(alpha = .3) +
                   geom_point(alpha = .6),
                 ggplot(stats, aes(x = S.obs, y = p_complete, ymin = p_complete_min, ymax = p_complete_max))+
                   labs(x = "Observed richness") +
                   scale_x_continuous(trans = 'sqrt') +
                   geom_errorbar(alpha = .3) +
                   geom_point(alpha = .6))
  
  plots <- c(lapply(plots1, function(p){
    p+
      labs(x = "Proportion of expected OTUs observed")+
      geom_vline(xintercept = threshold, col = "red")+
      theme_tufte()
  }),
  lapply(plots2, function(p){
    p+
      labs(y = "Proportion of expected OTUs observed")+
      geom_hline(yintercept = threshold, col = "red")+
      theme_tufte()
  }))
  
  print(plot_grid(plots[[1]]+
                    theme(axis.title.x = element_blank(),
                          axis.text.x = element_blank(),
                          axis.ticks.x = element_blank()), 
                  plots[[3]]+
                    theme(axis.title.x = element_blank(),
                          axis.text.x = element_blank(),
                          axis.ticks.x = element_blank()),
                  plots[[2]],
                  plots[[4]],
                  nrow = 2, align = "v")
  )
  
  output = list(
    minpass = row.names(stats[stats$p_complete_min >= threshold,]),
    meanpass = row.names(stats[stats$p_complete >= threshold,]),
    maxpass = row.names(stats[stats$p_complete_max >= threshold,])
  )
  
  return(list(passnames = output, stats = stats, plots = plots))
}

# Load data ---------------------------------------------------------------
  # Load the three core data tables and your tree, if you have one. 
  # If you have more than these data, great! Add loading commands for them
  # here. I would generally suggest merging and/or reconfiguring them as 
  # necessary to be similar to the configuration of these tables, or to
  # integrate their data within them.

  # Load the reads table output by vsearch --usearch_global --tabbedout

reads <- read.table("reads.tsv", header = T, sep = "\t", row.names = 1, comment.char = '')

    # Generally, samples should be rows and species should be columns, so
    # we transpose the usearch_global --tabbedout reads table

reads <- t(reads)

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
taxonomy <- read.table("SINTAXtaxonomy.tsv", sep = "\t", row.names = 1)
colnames(taxonomy) <- c("taxonomy", "strand", "selected")

    # From MEGAN, after using the script get_NCBI_taxonomy.py
taxonomy <- read.table("MEGANtaxonomy.tsv", sep = "\t", row.names = 1, header = T)

  # Load tree
tree <- read.tree("otus.nwk")

# Filtering low abundance OTU incidences ----------------------------------

# The filtering done in the metabarcoding pipeline removes whole ASVs that
# are potentially wrong, but doesn't do anything to filter putatively 
# incorrect incidences of valid OTUs from samples. This is an area we're 
# working on, but for now we do this by examining OTU incidence per sample

# Read number filtering 
  # Set counts less than a certain value to 0 - this isn't recommended
  # because it doesn't take into account total reads per sample. It also
  # generates difficulties with some of the accumulation tests later

#MBC_reads[MBC_reads < 2] <- 0

# Read proportion filtering
  # Set counts that are less than a certain proportion of reads in a 
  # sample to 0. 
  
  # First generate a filtering matrix, which gives the proportion of 
  # total reads in a sample for each OTU.
samplewise_p <- apply(reads, 1, function(s) s/sum(s) ) %>% t()

  # Could then explore this matrix. At the simplest level, produce a
  # histogram showing the distribution of proportion of per-OTU sample
  # incidences that exceed a threshold value (here 0.5%, i.e. 0.005).
apply(samplewise_p, 2, function(o) sum(o >= 0.005) / sum(o > 0)) %>% hist()

  # Could also find the maximum proportion for each OTU
maxprops <- apply(samplewise_p, 2, function(o) max(o)) %>% sort()

  # Can use this to show how many OTUs would be lost completely for 
  # different thresholds (because maxprops is sorted)
  # E.g. what value would retain 95% of OTUs:
maxprops[floor((1 - 0.95) * length(maxprops))]
  # E.g. what proportion of OTUs would be dropped if drop incidences less
  # than 0.5%
sum(maxprops <= 0.000125) / length(maxprops)

  # An important consideration is that no sample should be left without any
  # OTUs that occur only once, otherwise estimated richness cannot be 
  # calculated. Given that the minimum >0 proportion of any OTU within a 
  # sample is likely going to be for OTUs that only have 1 read, we can 
  # find the maximum threshold that fits this consideration
apply(samplewise_p, 1, function(s) min(s[s>0])) %>% min()

  # If this value is too low, you could pick a higher value and throw away
  # samples later. This isn't ideal though...

  # So, pick a threshold. This is very data-dependent. Think carefully
threshold <- 0.000125

  # Apply the threshold to the reads

reads[samplewise_p < threshold] <- 0

  # Crucially, we now clean up, dropping any samples that no longer have
  # any reads, and any OTUs that no longer have any reads. 
reads <- reads[rowSums(reads) > 0, colSums(reads) > 0]

rm(samplewise_p, maxprops, threshold)

# Organise data -----------------------------------------------------------
  # Here we make sure all the data corresponds properly to each other and 
  # any reconfiguration or merging is done as needed.

# READS

  # If we have multiple metabarcoding samples that actually represent the 
  # same ecological sample, here we might merge them together

## Some merging using rowsum(reads, group = XXX)

# METADATA

  # Next, check that our metadata and reads table correspond. For both 
  # tables, the sample names are in the row.names
    
    # Are all MBC samples in metadata? If not, problem!!!
row.names(reads) %in% row.names(metadata) 

    # Make the metadata correspond to the reads
metadata <- metadata[row.names(reads), ]
      # This keeps only metadata rows where the sample occurs in MBC, and 
      # sorts the rows to match the MBC row order

# TAXONOMY

  # The data in the taxonomy table needs to be filtered and separated out
    # As above, make sure all the OTUs are present then drop any taxonomy
    # records not present in the reads table
colnames(reads) %in% row.names(taxonomy)
taxonomy <- taxonomy[colnames(reads),]

  # Parse and reorganise a SINTAX taxonomy (skip this if using MEGAN)
{
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

    # Alternatively, We may only wish to include taxonomy info if it's above a certain score
#taxonomy <- dcast(taxdetailed[taxdetailed$score > 0.65, ], otu ~ level, value.var = "taxon", fill = '')
  
    # Set the row names
row.names(taxonomy) <- taxonomy$otu
taxonomy <- subset(taxonomy, select = -otu)    
}

  # Parse a MEGAN taxonomy
    # All we need to do here is drop the taxid column
taxonomy <- subset(taxonomy, select = -taxid)

    # Finally, re-sort the taxonomy
taxonomy <- taxonomy[colnames(reads), ]

rm(taxdetailed)

# TREE
  
  # Check that all of the OTUs are in the tree
all(row.names(reads) %in% tree$tip.label)

  # Drop any tips of the tree that are not in the OTUs - useful if you
  # used a scaffold or reference set
tree <- drop.tip(tree, tree$tip.label[! tree$tip.label %in% row.names(reads)])

  # Sort the OTUs according to the tree - will also drop any OTUs that are
  # not in the tree!
reads <- reads[tree$tip,]

  # Re-sort the taxonomy
taxonomy <- taxonomy[colnames(reads),]

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
reads <- reads[check$passnames$meanpass, ]

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


