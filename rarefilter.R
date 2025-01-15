# Setup -------------------------------------------------------------------

rm(list = ls())
options(stringsAsFactors = F)

# Load libraries ----------------------------------------------------------

suppressMessages(require(getopt))
suppressMessages(require(tidyr))
suppressMessages(require(dplyr))
suppressMessages(require(magrittr))
suppressMessages(require(purrr))
suppressMessages(require(ggplot2))
#suppressMessages(require(tibble))
#suppressMessages(require(knitr))


# Load functions ----------------------------------------------------------

plotter <- function(reads, norm, vline = NULL){
  caption <- "The lower half of the distribution of normalised read counts in your samples (bars), and the expected normal distribution (blue line).\n Potentially erroneous low counts will appear as bars above the normal line on the left side of the plot.\n"
  plot <- ggplot(mapping = aes(x = log10p)) + 
    geom_histogram(data = reads, binwidth = bw) + 
    geom_line(data = norm, aes(y = norm), col = "blue") + 
    scale_x_continuous(labels =  ~ signif(10^.x, digits = 2), breaks = scales::pretty_breaks(n = 10)) + 
    coord_cartesian(xlim = c(min(reads$log10p), max(vline, moments$meanlog10))) + 
    labs(x = "Normalised read counts", y = "Frequency") + 
    theme_bw()
  
  if(!is.null(vline)){
    plot + 
      geom_vline(xintercept = log10(vline), lty = "dashed", col = "red") + 
      labs(caption = paste(caption, "The vertical dashed red line represents your currently selected threshold, below which all counts will be dropped."))
  } else {
    plot + 
      labs(caption = paste(caption, "Based on this, select a threshold below which counts should be dropped. A strict option would be the x value\nat which the bars appear to meet the normal distribution; to be less conservative, choose a lower value."))
  }
}

# Set global variables ----------------------------------------------------



# Set up options ----------------------------------------------------------
# col1: long flag name
# col2: short flag name
# col3: 0 = no argument, 1 = required, 2 = optional
# col3: logical, integer, double, complex, character
# col5: optional, brief description

spec <- matrix(c(
  'help'           , 'h', 0, "logical"  , "show this helpful message.",
  'reads'          , 'r', 1, "character", "path to a tsv table recording the read count for each ASV within each sequenced sample, the first column (ASV_ID) must contain asv names and all subsequent columns samples.",
  'taxonomy'       , 'x', 2, "character", "path to a csv table with taxonomy data, with column header, from which to remove discarded ASVs. The first column will be assumed to be ASV IDs.",
  'mergetable'     , 'm', 2, "character", "optional path to a two-column csv that groups sequenced samples (column 1) into true samples (column 2), for example in the case of sequencing replicates. No column header.",
  'keepsingletons' , 's', 0, "logical"  , "turn off the default behaviour of removing all singletons no matter the threshold chosen",
  'outprefix'      , 'o', 2, "character", "prefix file path to write intermediate threshold plots and final filtered reads data"
), byrow = T, ncol = 5)

# Read options and do help -----------------------------------------------

opt <- getopt(spec)

if ( is.null(opt) | !is.null(opt$help) ){
  message(getopt(spec, usage = T))
  q(status = 1)
}

rm(spec)

opt$reads <- "ampliseqout_2024-02-20/dada2/ASV_table.tsv"

if( is.null(opt$keepsingletons) ){
  opt$keepsingletons <- F
}

# Load data and check ----------------------------------------------------

reads <- read.table(opt$reads, header = T)

if(names(reads)[1] != "ASV_ID"){
  stop(" the first column of the supplied reads tsv is not called \"ASV_ID\"")
}


# Report initial statistics -----------------------------------------------

ncols <- ncol(reads)-1
nrows <- nrow(reads)

message("Loaded initial reads data for ", nrows, " ASVs and ", ncols, 
        " sequenced samples: ", sum(as.matrix(reads[,-1])), " total reads.")

drop0cols <- sum(colSums(reads[,-1])==0)
drop0rows <- sum(rowSums(reads[,-1])==0)
message("Dropped ", drop0rows, " ASVs and ", drop0cols, " sequenced samples for having 0 reads.")
if(drop0cols > 0 | drop0rows > 0){
  message("Filtering ", nrows - drop0rows, "ASVs and ", ncols - drop0cols, " sequenced samples that have > 0 reads.")
}

# Reformat reads ----------------------------------------------------------

reads %<>% 
  rename_with(~ paste(1:ncols, .x), .cols = !ASV_ID) %>%
  pivot_longer(!ASV_ID, values_to = "count") %>%
  filter(count > 0) %>%
  separate(name, c("sampleid", "samplename"), sep = " ")
  
reformatsum <- reads %>% summarise_all(n_distinct)

if(ncols != reformatsum$samplename + drop0cols){
  warning("Some sequenced sample names are duplicated and these samples will be merged in a later step. The most common cause is sequencing the same sample on different runs with the same ID, in which case there is no cause for concern.")
}
if(nrows != reformatsum$ASV_ID + drop0rows){
  warning("Some ASV IDs are duplicated and these samples will be merged in a later step.")
}


# Initial plots -----------------------------------------------------------

# Distribution of ASVwise read counts
reads %>% 
  group_by(ASV_ID) %>%
  summarise(`Reads per ASV` = sum(count)) %>%
  ggplot(aes(x = `Reads per ASV`)) + 
  geom_histogram() + 
  scale_x_log10() + 
  theme_bw()

# Distribution of samplewise read counts
reads %>%
  group_by(sampleid) %>%
  summarise(`Reads per sequenced sample` = sum(count)) %>%
  ggplot(aes(x = `Reads per sequenced sample`)) + 
  geom_histogram() + 
  scale_x_log10() + 
  theme_bw()


# Normalised filter -------------------------------------------------------

message("Computing normalised read counts for filtering.")

# Calculate proportion of total reads in a sample for each ASV. As this is likely
# not normally distributed, log transform

reads %<>%
  group_by(sampleid) %>%
  mutate(p = count/sum(count),
         log10p = log10(p)) %>%
  ungroup()

# Calculate binwidth for plots
bw <- diff(range(reads$log10p))/60

# Find the moments of the normal distribution of log10p

moments <- reads %>%
  ungroup() %>%
  summarise(mean = mean(count),
            meanlog10 = mean(log10p),
            sdlog10 = sd(log10p))


# Construct an idealised distribution
norm <- 
  data.frame(log10p = seq(min(reads$log10p), 
                                max(reads$log10p), 
                                by = 0.01)) %>% 
  mutate(norm = dnorm(log10p, 
                      mean = moments$meanlog10, 
                      sd = moments$sdlog10),
         norm = norm/100 * bw * sum(reads$count),
         p = 10^log10p)


# Filter singletons -------------------------------------------------------

if(! opt$keepsingletons){
  reads <- reads %>%
    filter(count > 1)
  readsmultsum <- reads %>% summarise_all(n_distinct)
  message("Removing singleton records retains ", 
          readsmultsum$ASV_ID, "/", reformatsum$ASV_ID, " ASVs and ", 
          readsmultsum$sampleid, "/", reformatsum$sampleid, " sequenced samples.")
}

# Decide on filtering threshold and apply it ------------------------------

message("Please see rarefilter_thresholdchooser.pdf to review the normalised read counts. You may then experiment with the effect of different thresholds.")
ggsave(paste0(opt$outprefix, "rarefilter_thresholdchooser.pdf"), plot = plotter(reads, norm), device = "pdf", width = 8, height = 11)

choice <- 1
threshold <- 0

while(choice == 1){
  threshold <- readline("Please enter a threshold to test:")
  if( is.na(suppressWarnings(as.numeric(threshold))) ) {
    message("Error: threshold can't be interpreted as a number.")
    next
  } else {
    threshold <- as.numeric(threshold)
  }
  if(threshold < 0 | threshold > 1){
    message("Error: threshold must be between 0 and 1 inclusive.")
    next
  }
  filename <- paste0(opt$outprefix, "rarefilter_thresholdchooser_", threshold, ".pdf")
  ggsave(filename, 
         plot = plotter(reads, norm, threshold), 
         device = "pdf", width = 8, height = 11)
  threshreads <- reads %>% 
    filter(p > threshold)
  threshreadssum <- threshreads %>% summarise_all(n_distinct)
  message("Choosing a threshold of ", 
          threshold, " will retain ", 
          threshreadssum$ASV_ID, "/", reformatsum$ASV_ID, " ASVs and ", 
          threshreadssum$sampleid, "/", reformatsum$sampleid, " sequenced samples. ",
          "Please see ", filename, " to review this threshold on the histogram." )
  choice <- menu(c("Try another threshold", "Filter with this threshold"), title="Do you want to:")
}

message("Threshold filtering complete.")


# Merge samples -----------------------------------------------------------

mergedreads <- threshreads %>%
  group_by(ASV_ID, samplename) %>%
  summarise(count = sum(count)) %>%
  ungroup()

merge1sum <- mergedreads %>% summarise_all(n_distinct)

if(ncols != reformatsum$samplename + drop0cols | 
   nrows != reformatsum$ASV_ID + drop0rows) {
  message("Merging duplicate ASV ids and/or duplicate sequence sample names complete, ", merge1sum$samplename, " sequenced samples remain.")
}

if( !is.null(opt$mergetable) ){
  mergedetails <- read.csv(opt$mergetable, header = F) %>%
    setNames(c("sample", "samplename"))
  
  mergedreads <- left_join(mergedreads, mergedetails, by = c("samplename")) %>%
    group_by(ASV_ID, sample) %>%
    summarise(count = sum(count)) %>%
    ungroup() %>%
    rename("samplename" = "sample")
  
  merge2sum <- mergedreads %>% summarise_all(n_distinct)
  
  message("Merging using the supplied mergetable complete, ", merge2sum$samplename, " samples remain.")

}

# Convert and output reads ------------------------------------------------

mergedreads %>%
  pivot_wider(names_from = samplename, values_from = count, values_fill = 0) %>%
  write.csv(paste0(opt$outprefix, "filtered_readstable.csv"), row.names = F, quote = F)


# Filter and reoutput taxonomy --------------------------------------------

if( !is.null(opt$taxonomy) ){
  read.csv(opt$taxonomy) %>%
    filter(.[[1]] %in% mergedreads$ASV_ID) %>%
    write.csv(paste0(opt$outprefix, "filtered_taxonomy.csv"), row.names = F, quote = F)
}

