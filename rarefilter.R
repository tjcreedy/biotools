# Setup -------------------------------------------------------------------

rm(list = ls())
options(stringsAsFactors = F)

# Load libraries ----------------------------------------------------------

suppressMessages(require(getopt))
suppressMessages(require(tidyr))
suppressMessages(require(dplyr))
suppressMessages(require(magrittr))
suppressMessages(require(ggplot2))


# Load functions ----------------------------------------------------------

p <- paste0
vline = NULL
trim = TRUE
plotter <- function(reads, vline = NULL, trim = TRUE){
  caption <- p("The lower half of the distribution of normalised read counts in your samples (bars), and the expected ",
               "normal distribution (blue line).\n Potentially erroneous low counts will appear as bars above the ",
               "normal line on the left side of the plot.\n")
  
  moments <- reads %>%
    ungroup() %>%
    summarise(mean = mean(count),
              meanlog10 = mean(log10p),
              sdlog10 = sd(log10p))
  
  plot <- ggplot(reads, aes(x = log10p)) + 
    geom_histogram(aes(y = after_stat(density)), bins = 60) + 
    stat_function(fun = dnorm, args = list(mean = moments$meanlog10, sd = moments$sdlog10), col = "blue") + 
    scale_x_continuous(labels =  ~ signif(10^.x, digits = 2) %>% scales::comma(), breaks = scales::pretty_breaks(n = 10)) + 
    labs(x = "Normalised read counts", y = "Density") + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  
  if(trim){
    xlimmax <- max(if(is.null(vline)) -Inf else log10(vline), mean(reads$log10p))
    plot <- plot + coord_cartesian(xlim = c(min(reads$log10p), xlimmax))
  }
  
  if(!is.null(vline)){
    plot + 
      geom_vline(xintercept = log10(vline), lty = "dashed", col = "red") + 
      labs(caption = paste(caption, "The vertical dashed red line represents your currently selected threshold, below",
                           "which all counts will be dropped."))
  } else {
    plot + 
      labs(caption = paste(caption, "Based on this, select a threshold below which counts should be dropped. A strict",
                           "option would be the x value\nat which the bars appear to meet the normal distribution; to",
                           "be less conservative, choose a lower value."))
  }
}

typeline <- function(msg="Enter text: ") {
  if (interactive() ) {
    txt <- readline(msg)
  } else {
    cat(msg);
    txt <- readLines("stdin",n=1);
  }
  return(txt)
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
  'reads'          , 'r', 1, "character", p("path to a tsv table recording the read count for each ASV within each ",
                                            "sequenced sample, the first column (ASV_ID) must contain asv names and ",
                                            "and all subsequent columns samples."),
  'taxonomy'       , 'x', 2, "character", p("path to a csv table with taxonomy data, with column header, from which to", 
                                            " remove discarded ASVs. The first column will be assumed to be ASV IDs."),
  'mergetable'     , 'm', 2, "character", p("optional path to a two-column csv that groups sequenced samples (column ",
                                            "1) into true samples (column 2), for example in the case of sequencing ", 
                                            "replicates. No column header."),
  'threshold'      , 't', 2, "double"   , p("optionally supply a threshold to filter without reviewing graphs (not ",
                                            "recommended)"),
  'keepsingletons' , 's', 0, "logical"  , p("turn off the default behaviour of removing all singletons no matter the ",
                                            "threshold chosen"),
  'outprefix'      , 'o', 2, "character", p("prefix file path to write intermediate threshold plots and final ",
                                            "filtered reads data")
), byrow = T, ncol = 5)

# Read options and do help -----------------------------------------------

opt <- getopt(spec)

if ( is.null(opt) | !is.null(opt$help) ){
  message(getopt(spec, usage = T))
  q(status = 1)
}

rm(spec)

#opt$reads <- "ASV_codon_filtered.table.tsv"
#opt$taxonomy <- "b2t_out_taxonomy.csv"
#opt$threshold <- 3e-05

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

message("\nLoaded initial reads data for ", nrows, " ASVs and ", ncols, 
        " sequenced samples: ", sum(as.matrix(reads[,-1])), " total reads.")

drop0cols <- sum(colSums(reads[,-1])==0)
drop0rows <- sum(rowSums(reads[,-1])==0)
message("Dropped ", drop0rows, " ASVs and ", drop0cols, " sequenced samples for having 0 reads.")
if(drop0cols > 0 | drop0rows > 0){
  message("Filtering ", nrows - drop0rows, " ASVs and ", ncols - drop0cols, " sequenced samples that have >0 reads.")
}

# Reformat reads ----------------------------------------------------------

reads %<>% 
  rename_with(~ paste(1:ncols, .x), .cols = !ASV_ID) %>%
  pivot_longer(!ASV_ID, values_to = "count") %>%
  filter(count > 0) %>%
  separate(name, c("sampleid", "samplename"), sep = " ")
  
reformatsum <- reads %>% summarise_all(n_distinct)

if(ncols != reformatsum$samplename + drop0cols){
  warning("\nSome sequenced sample names are duplicated and these samples will be merged in a later step. The most ", 
          "common cause is sequencing the same sample on different runs with the same ID, in which case there is no ",
          "cause for concern.")
}
if(nrows != reformatsum$ASV_ID + drop0rows){
  warning("\nSome ASV IDs are duplicated and these samples will be merged in a later step.")
}


# Initial plots -----------------------------------------------------------

# Distribution of ASVwise read counts
# reads %>% 
#   group_by(ASV_ID) %>%
#   summarise(`Reads per ASV` = sum(count)) %>%
#   ggplot(aes(x = `Reads per ASV`)) + 
#   geom_histogram() + 
#   scale_x_log10() + 
#   theme_bw()
# 
# # Distribution of samplewise read counts
# reads %>%
#   group_by(sampleid) %>%
#   summarise(`Reads per sequenced sample` = sum(count)) %>%
#   ggplot(aes(x = `Reads per sequenced sample`)) + 
#   geom_histogram() + 
#   scale_x_log10() + 
#   theme_bw()


# Filter singletons -------------------------------------------------------

if(! opt$keepsingletons){
  reads <- reads %>%
    filter(count > 1)
  readsmultsum <- reads %>% summarise_all(n_distinct)
  message("Removed singleton records, ", 
          readsmultsum$ASV_ID, "/", reformatsum$ASV_ID, " ASVs and ", 
          readsmultsum$sampleid, "/", reformatsum$sampleid, " sequenced samples remain.")
}

# Normalised filter -------------------------------------------------------

message("\n\nComputing normalised read counts for filtering.")

# Calculate proportion of total reads in a sample for each ASV. As this is likely
# not normally distributed, log transform

reads %<>%
  group_by(sampleid) %>%
  mutate(p = count/sum(count),
         log10p = log10(p)) %>%
  ungroup()


# Decide on filtering threshold and apply it ------------------------------

if(! is.null(opt$threshold) ){
  threshreads <- reads %>% 
    filter(p > opt$threshold)
  threshreadssum <- threshreads %>% summarise_all(n_distinct)
  message("\nFiltering at a threshold of ", 
          opt$threshold, " retained ", 
          threshreadssum$ASV_ID, "/", reformatsum$ASV_ID, " ASVs and ", 
          threshreadssum$sampleid, "/", reformatsum$sampleid, " sequenced samples.")
} else {
  message("\nPlease see rarefilter_thresholdchooser.pdf to review the normalised read counts. You may then experiment ",
          "with the effect of different thresholds.")
  ggsave(paste0(opt$outprefix, "rarefilter_thresholdchooser.pdf"), 
         plot = plotter(reads), device = "pdf", width = 8, height = 11)
  
  keeptrying <- TRUE
  threshold <- 0
  
  while(keeptrying){
    threshold <- typeline("\nPlease enter a threshold to test: ")
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
           plot = plotter(reads, threshold), 
           device = "pdf", width = 8, height = 11)
    threshreads <- reads %>% 
      filter(p > threshold)
    threshreadssum <- threshreads %>% summarise_all(n_distinct)
    message("\nChoosing a threshold of ", 
            threshold, " will retain ", 
            threshreadssum$ASV_ID, "/", reformatsum$ASV_ID, " ASVs and ", 
            threshreadssum$sampleid, "/", reformatsum$sampleid, " sequenced samples. ",
            "Please see ", filename, " to review this threshold on the histogram." )
    while(1){
      choice <- typeline("\nDo you want to (T)ry another threshold or (F)ilter with this threshold? ")
      if(choice %in% c("F", "f", "T", "t")){
        break
      } else {
        message("Input not understood.")
      }
    }
    keeptrying <- as.logical(toupper(choice))
  }
}

message("\nThreshold filtering complete, finalising outputs.")

# Merge samples -----------------------------------------------------------

mergedreads <- threshreads %>%
  group_by(ASV_ID, samplename) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  ungroup()

merge1sum <- mergedreads %>% summarise_all(n_distinct)

if(ncols != reformatsum$samplename + drop0cols | 
   nrows != reformatsum$ASV_ID + drop0rows) {
  message("\nMerging duplicate ASV ids and/or duplicate sequence sample names complete, ", merge1sum$samplename, 
          " sequenced samples remain.")
}

if( !is.null(opt$mergetable) ){
  mergedetails <- read.csv(opt$mergetable, header = F) %>%
    setNames(c("sample", "samplename"))
  
  mergedreads <- left_join(mergedreads, mergedetails, by = c("samplename")) %>%
    group_by(ASV_ID, sample) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    ungroup() %>%
    rename("samplename" = "sample")
  
  merge2sum <- mergedreads %>% summarise_all(n_distinct)
  
  message("\nMerging using the supplied mergetable complete, ", merge2sum$samplename, " samples remain.")

}

# Convert and output reads ------------------------------------------------

mergedreads %>%
  pivot_wider(names_from = samplename, values_from = count, values_fill = 0) %>%
  write.csv(paste0(opt$outprefix, "filtered_readstable.csv"), row.names = F, quote = F)


# Filter and reoutput taxonomy --------------------------------------------

if( !is.null(opt$taxonomy) ){
  read.csv(opt$taxonomy) %>% as_tibble()
    filter(.[[1]] %in% mergedreads$ASV_ID) %>%
    write.csv(paste0(opt$outprefix, "filtered_taxonomy.csv"), row.names = F, quote = F)
}

message("\nAll done!")
