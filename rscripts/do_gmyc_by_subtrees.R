#!/usr/bin/env Rscript

# Set options -------------------------------------------------------------

options(stringsAsFactors = F)

# Load libraries ----------------------------------------------------------

suppressMessages(require(getopt))
suppressMessages(require(ape))
suppressMessages(require(splits))
suppressMessages(require(foreach))
suppressMessages(require(doMC))

# Define functions --------------------------------------------------------

listdescendants <- function(tree, n, nodes = T, tips = T, inc.n = F){
  if(nodes == F & tips == F){
    stop("Nothing to return!")
  }
  if(n < 1 | n > Ntip(tree) + tree$Nnode){
    stop("Argument to n is not a valid node or tip")
  }
  
  chs <- tree$edge[tree$edge[,1] == n, 2]
  
  out <- unlist(sapply(chs, function(ch){
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

get_subtree_nodes <- function(tree, n, node = Ntip(tree) + 1){
  if(node < 0 | node > Ntip(tree) + tree$Nnode){
    stop("Argument to node is not a valid node")
  }
  if(node <= Ntip(tree)){
    return(node)
  }
  ndesc <- length(listdescendants(tree, n = node, nodes = F, tips = T))
  if(ndesc <= n){
    return(node)
  } else {
    chs <- tree$edge[tree$edge[,1] == node, 2]
    out <- sapply(chs, function(ch) get_subtree_nodes(tree = tree, n = n, node = ch))
    return(unlist(out))
  }
}

get_subtrees <- function(tree, n = NULL, k = NULL){
  if(is.null(n) & is.null(k)){
    stop("Supply a value to n or k")
  }
  if(! is.null(n) & ! is.null(k)){
    stop("Supply only n or k")
  }
  if(is.null(n)){
    n = floor(Ntip(tree)/k)
  }
  nodes = get_subtree_nodes(tree, n)
  subtrees <- lapply(nodes, function(node) {
    keep.tip(tree, listdescendants(tree, node, nodes = F, inc.n = T))
  })
  return(subtrees)
}

# Set up and read options -------------------------------------------------

spec <- matrix(c(
  'tree'     , 'n', 1, 'character',
  'outdir'   , 'o', 1, 'character',
  'nsubtrees', 'k', 2, 'integer',
  'subsize'  , 's', 2, 'integer',
  'threads'  , 't', 2, 'integer',
  'help'     , 'h', 0, 'logical'
), byrow = T, ncol = 4)


opt <- getopt(spec)

opt$tree = "upgma_tree.nwk"
opt$threads = 2
opt$outdir = "test/"
opt$nsubtrees = 4

# Parse options -----------------------------------------------------------

if ( !is.null(opt$help) ){
  cat(getopt(spec, usage = T))
  q(status = 1)
}

if ( !is.null(opt$nsubtrees) & !is.null(opt$subsize) ){
  stop("Give only one of nsubtrees or subsize")
}
if ( is.null(opt$nsubtrees) & is.null(opt$subsize) ){
  stop("Give one of nsubtrees or subsize")
}
if ( is.null(opt$threads) ){
  opt$threads <- 1
}

registerDoMC(opt$threads)

# Load data ---------------------------------------------------------------

tree <- read.tree(opt$tree)

# Get subtree list --------------------------------------------------------

subtrees <- list()
if( !is.null(opt$nsubtrees) ){
  subtrees <- get_subtrees(tree, k = opt$nsubtrees)
} else {
  subtrees <- get_subtrees(tree, n = opt$subsize)
}

# Extract subtrees with only 1 tip ----------------------------------------

subtrees1 <- list()
subtreesrun <- list()
for(tr in subtrees){
  if(tr$Nnode == 1){
    subtrees1[[length(subtrees1)+1]] <- tr
  } else {
    subtreesrun[[length(subtreesrun)+1]] <- tr
  }
}

# Perform GMYC on subtrees in parallel ------------------------------------

gmycresults <- foreach(i = 1:length(subtreesrun)) %dopar% gmyc(subtreesrun[[i]])

saveRDS(gmycresults, file = paste0(opt$outdir, "/allgmycresults.RDS"))

# Concatenate and output statistics and groupings -------------------------

summaries <- do.call("rbind", lapply(1:length(gmycresults), function(i) {
  sumgmyc <- capture.output(summary(gmycresults[[i]]))
  sumdat <- do.call("rbind", strsplit(sumgmyc, "\t", fixed = T))[-1, -1]
  names <- c("subtree", gsub(':', '', sumdat[,1]))
  dat <- c(i, sumdat[,2])
  return(matrix(dat, nrow = 1, dimnames = list(NULL, names)))
}))
write.csv(summaries, paste0(opt$outdir, "/summaries.csv"), row.names = F, quote = F)
rm(summaries)

speclist <- do.call("rbind", lapply(1:length(gmycresults), function(i) {
  cbind(subtree = i, spec.list(gmycresults[[i]]))
}))
speclist1 <- do.call("rbind", lapply(1:length(subtrees1), function(i){
  c(subtree = length(subtreesrun) + 1,
    GMYC_spec = 1,
    sample_name = subtrees1[[i]]$tip.label)
}))

write.csv(rbind(speclist, speclist1), file = paste0(opt$outdir, "/speclist.csv"), row.names = F, quote = F)
