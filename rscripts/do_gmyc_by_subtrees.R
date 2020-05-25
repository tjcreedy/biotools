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
  return(do.call(c, subtrees))
}

# Set up and read options -------------------------------------------------

spec <- matrix(c(
  'tree'     , 'n', 1, 'character',
  'outdir'   , 'o', 1, 'character',
  'subsize'  , 's', 2, 'integer',
  'threads'  , 't', 2, 'integer',
  'help'     , 'h', 0, 'logical'
), byrow = T, ncol = 4)


opt <- getopt(spec)

# opt$tree = "~/Documents/NHM_postdoc/Supervision/ZichenZ/demo_upgma_zz.nwk"
# opt$threads = 2
# opt$outdir = "./"
# opt$subsize = 25

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

message("Starting subtree division")

subtrees <- NULL
if( !is.null(opt$nsubtrees) ){
  subtrees <- get_subtrees(tree, k = opt$nsubtrees)
} else {
  subtrees <- get_subtrees(tree, n = opt$subsize)
}

message(paste0("Generated ", length(subtrees), " subtrees"))

# Extract subtrees with only 1 tip ----------------------------------------

onenode <- sapply(subtrees, function(tr) tr$Nnode == 1)
subtrees1 <- subtrees[onenode]
subtreesrun <- subtrees[!onenode]

write.tree(subtrees1, paste0(opt$outdir, "/subtrees_singletons.nwk"))
write.tree(subtreesrun, paste0(opt$outdir, "/subtrees_forgmyc.nwk"))

# Perform GMYC on subtrees in parallel ------------------------------------

message(paste0("Running GMYC on ", length(subtreesrun), " subtrees with >2 tips"))

gmycresults <- foreach(i = 1:length(subtreesrun)) %dopar% tryCatch(gmyc(subtreesrun[[i]]), error = function(x) NULL)
gmycresults <- gmycresults[! sapply(gmycresults, is.null)]

message(paste0("GMYC succeeded for ", length(gmycresults), " subtrees"))

saveRDS(gmycresults, file = paste0(opt$outdir, "/allgmycresults.RDS"))

# Concatenate and output statistics and groupings -------------------------
speclist <- do.call("rbind", lapply(1:length(gmycresults), function(i) {
  tryCatch(cbind(subtree = i, spec.list(gmycresults[[i]])), 
           error = function(x) NULL)
}))

slsucceed <- sum(1:length(gmycresults) %in% speclist$subtree)
message(paste0("Species delimitation succeeded for ", slsucceed, " GMYC results"))
speclist1 <- do.call("rbind", lapply(1:length(subtrees1), function(i){
  tips <- subtrees1[[i]]$tip.label
  data.frame(subtree = length(subtreesrun) + i,
             GMYC_spec = 1:length(tips),
             sample_name = tips)
}))
write.csv(rbind(speclist, speclist1), file = paste0(opt$outdir, "/speclist.csv"), row.names = F, quote = F)

summaries <- do.call("rbind", lapply(1:length(gmycresults), function(i) {
  sumgmyc <- capture.output(summary(gmycresults[[i]]))
  sumdat <- do.call("rbind", strsplit(sumgmyc, "\t", fixed = T))[-1, -1]
  names <- c("subtree", gsub(':', '', sumdat[,1]))
  dat <- c(i, sumdat[,2])
  return(matrix(dat, nrow = 1, dimnames = list(NULL, names)))
}))
write.csv(summaries, paste0(opt$outdir, "/summaries.csv"), row.names = F, quote = F)
rm(summaries)