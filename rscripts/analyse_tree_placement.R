#!/usr/bin/env Rscript

# Setup -------------------------------------------------------------------

rm(list = ls())
options(stringsAsFactors = F)

# Load libraries ----------------------------------------------------------

suppressMessages(require(getopt))
suppressMessages(require(geiger))
suppressMessages(require(ape))
suppressMessages(require(plyr))
suppressMessages(require(vegan))
suppressMessages(require(ggplot2))



# Load functions ------------------------------------------------------------------------------

listdescendants <- function(tree, n, nodes = T, tips = T, inc.n = F){
  if(nodes == F & tips == F){
    stop("Nothing to return!")
  }
  if(n < 1 | n > Ntip(tree) + tree$Nnode){
    stop("Argument to n is not a valid node or tip")
  }
  # if(inc.n == T){
  #   if(nodes == F & ! n <= Ntip(tree)) warning("inc.n = T ignored because n is an internal node and nodes = F")
  # }
  
  chs <- tree$edge[tree$edge[,1] == n, 2]
  
  out <- unlist(lapply(chs, function(ch){
    if(ch <= Ntip(tree)){
      if(tips == T){
        return(ch)
      } else if(tips == "labels"){
        return(tree$tip.label[ch])
      }
    } else {
      ot <- listdescendants(tree, ch, nodes, tips)
      if(nodes == T){
        ot <- c(ch,ot)
      }
      return(ot)
    }
  }))
  
  if(inc.n == T & (nodes == T | n <= Ntip(tree))){
    out <- c(n, out)
  }
  return(out)
}

find_monophyletic_subtrees <- function(tree, tips, start = Ntip(tree)+1){
  for(tip in tips){
    if(!tip %in% tree$tip.label){
      stop(paste("Tip", tip, "not found in the supplied tree"))
    }
  }
  
  currtips <- NULL
  if(start %in% 1:Ntip(tree)){
    currtips <- start
  } else {
    currtips <- listdescendants(tree = tree, n  = start, nodes = F, tip = T, inc.n = F)
  }
  
  currtips
  
  intips <- tree$tip.label[currtips] %in% tips
  
  if(all(intips)){
    return(c(start))
  } else if(! any(intips)){
    return(NULL)
  } else if(sum(intips) == 1){
    return(currtips[intips])
  } else {
    chs <- tree$edge[tree$edge[,1] == start, 2]
    return(unlist(sapply(chs, function(ch) find_monophyletic_subtrees(tree, tips, start = ch))))
  }
}

find_largest_outgroup_parent <- function(tree, tips){
  subtreenodes <- find_monophyletic_subtrees(tree, tips)
  subtreelength <- lapply(subtreenodes, function(n){
    length(listdescendants(tree, n, nodes = F))
  })
  return(subtreenodes[which.max(subtreelength)])
}

root_outgroup_fuzzy <- function(tree, outgroup){
  return(ladderize(root(tree, node = find_largest_outgroup_parent(tree, outgroup), resolve.root = T)))
}


sortgenes <- function(g){
  if( length(g) == length(geneorder)){
    return(geneorder)
  }
  if( length(g) == 1){
    return(g)
  }
  # Find sets of matches
  io <- match(c(geneorder, geneorder), g)
  ioNonmissingi <- which(!is.na(io))
  sets <- lapply(1:(length(ioNonmissingi)-1), function(a){
    i <- ioNonmissingi[a]
    later <- ioNonmissingi[(a+1):length(ioNonmissingi)]
    later <- later[(later-i) < 13]
    lapply(later, function(l) io[i:l])
  })
  sets <- do.call(c, sets)
  # Filter out sets that don't contain all input genes
  sets <- sets[sapply(sets, function(set) all(unique(na.omit(io)) %in% set))]
  # Find lengths
  setlengths <- sapply(sets, length)
  # Return genes corresponding to first set that is the shortest
  return(g[sets[setlengths == min(setlengths)][[1]]])
}

consec_genes <- function(a, b) diff(match(c(a, b), geneorder)) == 1
consec_genes <- Vectorize(consec_genes, USE.NAMES = F)


# Global variables ----------------------------------------------------------------------------

genes <- c("ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6")
geneorder <- c("ND2", "COX1", "COX2", "ATP8", "ATP6", "COX3", "ND3", "ND5", "ND4", "ND4L", "ND6", "CYTB", "ND1") 


# Set up options ----------------------------------------------------------
# col1: long flag name
# col2: short flag name
# col3: 0 = no argument, 1 = required, 2 = optional
# col3: logical, integer, double, complex, character
# col5: optional, brief description

spec <- matrix(c(
  'help'        , 'h', 0, "logical"  , "show this helpful message",
  'genetrees'   , 't', 1, "character", "path to a directory containing gene trees",
  'reference'   , 'r', 1, "character", "path to a file detailing the reference tips in the gene 
                                        trees, can either be a text file with tip names one per 
                                        line, or a phylogeny from which tip names will be read",
  'groups'      , 'g', 1, "character", "path to a file detailing the groups of tips to compare, each
                                        line should be a separated list of tips separated by tab, 
                                        space or comma",
  'outgroup'    , 'o', 1, "character", "path to a file specifying the outgroup tip(s) using which 
                                        the tree will be re-rooted"
), byrow = T, ncol = 5)
spec[,5] <- gsub("\n[\t ]*", "", spec[,5])

# Read options and do help -----------------------------------------------

opt <- getopt(spec)

setwd("/home/thomas/work/iBioGen_postdoc/MMGdatabase/phylogeny/duplicatechecking/")
opt$genetrees <- "5_trees"
opt$reference <- "4_backbone.tre"
opt$groups <- "31_duplicate_groups.txt"
opt$outgroup <- "outgroup"

if ( is.null(opt) | !is.null(opt$help) ){
  message(getopt(spec, usage = T))
  q(status = 1)
}

rm(spec)


# Parse trees ---------------------------------------------------------------------------------

# Find tree files and check numbers
treefiles <- list.files(opt$genetrees, full.names = T)
genefiles <- lapply(genes, function(g) {
  treefiles[grepl(paste0("[^A-Za-z0-9]", g, "[^A-Za-z0-9]"), treefiles)]
})
genefiles <- setNames(genefiles, genes)
nmatches <- sapply(genefiles, length)
if( any(nmatches > 1) ){
  stop("Error: found more than one tree for ", paste0(genes[nmatches > 1], collapse = ","))  
}

# Read trees
genetrees <- lapply(genefiles[nmatches == 1], read.tree)

# Read outgroup
outgroup <- readLines(opt$outgroup)

# Check outgroup
outgroupmissing <- sapply(genetrees, function(tr) all(!outgroup %in% tr$tip.label))
if( any(outgroupmissing) ){
  stop("Error: all outgroup terminal(s) missing from some gene trees: ", 
       paste(names(genetrees)[outgroupmissing], collapse = ","))
}

# Reroot
genetrees <- lapply(genetrees, function(tr) root_outgroup_fuzzy(tr, outgroup = outgroup[outgroup %in% tr$tip.label]))


# Parse reference list ------------------------------------------------------------------------

# Read in references from text or tree
ref <- readLines(opt$reference)
if(length(ref) == 1){
  ref <- read.tree(opt$reference)$tip.label
}

# Check which references are missing and drop those not present in all trees
checkpresence <- lapply(genetrees, function(gt){
  ref[! ref %in% gt$tip.label]
})
absent <- unique(unlist(checkpresence))
ref <- ref[ !ref %in% absent ]

# Report to use
if( length(absent) > 0 ){
  msg <- paste("some references not present in one or more gene trees: \n", 
               paste(absent, collapse = ', '), 
               "\n", length(ref), "references remain.")
  if( length(ref) > 1 ){
    message("Warning: ", msg)
  } else {
    stop("Error: ", msg)
  }
}

# Parse groups --------------------------------------------------------------------------------

groups <- strsplit(readLines(opt$groups), "[\t ,]")


# Drop unneeded tips from the trees -----------------------------------------------------------

reqtips <- c(unique(unlist(groups)), ref)

genetrees <- lapply(genetrees, function(gt) keep.tip(gt, reqtips[reqtips %in% gt$tip.label]))


# Compute and compile distances ---------------------------------------------------------------

subset_trees <- function(trees, focctg, ref, branchlengths = 0.1){
  # Remove trees for genes not present
  out <- trees[sapply(trees, function(tree) any(focctg %in% tree$tip.label))]
  # Drop other tips
  out <- lapply(out, function(ot) drop.tip(ot, ot$tip.label[!ot$tip.label %in% c(ref, focctg)]))
  # Drop branch lengths
  if( branchlengths ){
    out <- lapply(out, function(ot) {
      ot$edge.length <- rep(branchlengths, length(ot$edge.length))
      return(ot)
    })
  }
  return(out)
}

genedata <- lapply(unlist(groups), function(focctg){
  
  # Subset trees and compute distances
  placetrees <- subset_trees(genetrees, focctg, ref)
    # Compute and compile distances
  genedists <- lapply(placetrees, cophenetic.phylo)
  genedists <- lapply(names(genedists), function(gene){
    data.frame(contig = focctg, gene = gene, genedists[[gene]][focctg, ref, drop = F])
  })
  genedists <- do.call(rbind, genedists)
  row.names(genedists) <- NULL
  
  return(genedists)
})
genedata <- do.call(rbind, genedata)


# Sort genes in each contig
genedata <- ddply(genedata, ~contig, function(gdc){
  return(gdc[match(na.omit(sortgenes(gdc$gene)), gdc$gene), ])
})


# Do ordinations ------------------------------------------------------------------------------

# Ordinate and extract scores
genedists <- genedata[, ref]
pca <- princomp(genedists)
genepoint <- cbind(genedata[, c("contig", "gene")], pca$scores)

# Prep for plotting
axes <- 1:2
pnts <- genepoint[genepoint$contig == genepoint$contig[1], ]
geneplot <- ddply(genepoint, ~ contig, function(pnts){
  # Make dataframe
  dummy <- matrix(nrow = 1, ncol = 3, dimnames = list(NULL, names(pnts[2:4])))
  if(nrow(pnts) > 1){
    out <- data.frame(pnts[, c(1:2, 2+axes)],
                      rbind(pnts[2:nrow(pnts), c(2, 2+axes)], dummy)
    )
  } else {
    out <- cbind(pnts[, c(1:2, 2+axes)], dummy)
  }
  names(out)[2:7] <- paste0(rep(c("", "next"), each = 3), c("gene", "x", "y"))
  
  # Calculate consecutive and distance
  if( nrow(pnts) == 1){
    return(data.frame(out, consecutive = NA, distance = NA))
  } else {
    out$consecutive <- consec_genes(out$gene, out$nextgene)
    out$distance <- c(mapply(function(i, j) dist(pnts[c(i,j), -c(1:2)]), 1:(nrow(out)-1), 2:nrow(out)), NA)
    return(out)
  }
})

geneplot$genepair <- factor(paste(geneplot$gene, geneplot$nextgene, sep = "-"), 
                            levels = paste(geneorder, geneorder[c(2:13, 1)], sep = "-"))



# Create group reports ------------------------------------------------------------------------

roundaway <- function(x) sign(x)*ceiling(abs(x))
palette <- c('#1f78b4','#33a02c','#e31a1c','#ff7f00','#6a3d9a','#b15928',
             '#a6cee3','#b2df8a','#fb9a99','#fdbf6f','#cab2d6','#ffff99')


library(cowplot)


plots <- lapply(groups, function(grpmem){
  
  # Set up subset for plotting
  contigplot <- geneplot[geneplot$contig %in% grpmem, ]
  # Create nudging vector
  contigplot$nudge <- as.numeric(as.factor(contigplot$contig))
  contigplot$nudge <- contigplot$nudge-mean(contigplot$nudge)
  
  # Create palette
  cols <- setNames(palette[seq_along(grpmem)], grpmem)
  
  vnudgeweight = .2
  vplot <- ggplot(geneplot[!is.na(geneplot$nextgene) & geneplot$consecutive, ], 
                  aes(x = genepair, y = distance)) + 
    geom_violin() +
    geom_point(data = na.omit(contigplot), aes(col = contig, y = distance + nudge * vnudgeweight)) + 
    geom_line(data = na.omit(contigplot), aes(col = contig, group = contig, y = distance + nudge * vnudgeweight)) + 
    scale_colour_manual(guide = "none", values = cols) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  onudgeweight = diff(range(geneplot$x))/250
  oplot <- ggplot(contigplot, aes(col = contig)) + 
    geom_text(aes(x = x + nudge * onudgeweight, y = y, label = gene)) + 
    geom_segment(aes(x = x + nudge * onudgeweight, y = y, 
                     xend = nextx + nudge * onudgeweight, yend = nexty, 
                     linetype = consecutive), 
                 alpha = .7)+
    scale_colour_manual(values = cols) +
    scale_linetype(guide = "none") + 
    lims(x = 10*roundaway(range(geneplot$x)/10),
         y = 10*roundaway(range(geneplot$y)/10)) + 
    theme_bw() + 
    theme(legend.position = c(.5, .85), axis.title = element_blank())
  
  
  contigtrees <- subset_trees(genetrees, grpmem, ref)
  par(mfrow = c(3, 5), mar = c(0,0,0,0))
  invisible(lapply(geneorder, function(gene){
    if(gene %in% names(contigtrees)){
      gt <- contigtrees[[gene]]
      plot.phylo(gt, show.tip.label = F)
      title(gene, line = -8, adj = .2)
      tiplabels(grpmem, match(grpmem, gt$tip.label), frame = "none", adj = 0, offset = .05, 
                col = cols)
    } else plot.new()
  }))
  invisible(lapply(1:2, function(i) plot.new()))
  tplot <- recordPlot()
  outplot <- plot_grid(tplot, 
                       plot_grid(oplot, vplot, ncol = 2), 
                       nrow = 2, rel_heights = c(3, 2))
  # pdf("test.pdf", width = 12, height = 12)
  # outplot
  # dev.off()
  return(outplot)
})

pdf("placement_report.pdf", width = 12, height = 12)
par(mfrow = c(1,1), mar = c(0,0,0,0))
plot.new()
text(x = .025, y = .5, 
     "Tree placement report


The following plots report the placement of contigs relative to reference sequence for each gene.

(Top panel:) A tree is shown for each gene represented by the contigs under comparison, displaying 
only the terminals for the contigs and the reference set (all other terminals removed).

(Bottom left panel:) The ordination illustrates the position of each gene on each contig relative 
to the reference set. The closer genes are, the more consistent their placement. Lines join 
consecutive genes from the same contig.

(Bottom right panel:) The violin plot reports the \"distance\" (i.e. dissimilarity) between the
placement of consecutive genes on the same contig (i.e. equivalent to the length of the lines 
between genes on the ordination). The violin plots report the distribution from the input dataset, 
points report the values for the contigs under comparison, lines join contigs. Low values mean that
both genes were placed in similar locations on their repsective genes, high values the opposite. 
High values may indicate a chimeric contig joined at some point between or within one of the 
consecutive genes.

NB: points and lines on bottom panels are printed with slight offsets to better show colliding 
data.
", adj = 0, cex = .9)
plots
dev.off()

# This basic pipeline could be used to detect chimeras automatically
# Inputs: reference set of known good, and contigs to consider
# Process: 
# 1. align all quickly and look for disagreements on circular contig location (e.g. two 
# circularised contigs with different origins)
# 2. break up contigs that disagree to ensure harmonious alignment
# 3. align again more accurately (possibly create reference alignment/consensus and then break up
# focal sequence to kmers and map each kmer)
# 4. use sliding window across alignment to assess similarity of focal sequences to references (
# or look at individual kmer results)
# 5. Use the similarity/dissimilarity against references in place of tree placement similarity
# in above logic - look at pattern of placement similarity across contig to find breaks