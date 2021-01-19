# prep for notung reconciliation and rooting
# thereafter run grampa

library("ape")
library("dplyr")
library("parallel")

rm(list=ls())

# function to remove _, @, and flip the names to postorder
renametips <- . %>% gsub("_",".",.) %>% strsplit(., "@") %>% do.call(rbind, .) %>% data.frame(., stringsAsFactors=F) %>% mutate(X3=paste(X2,"_",X1,sep=""))

# function to remove outgroup pelagophytes
dropOutgroup <- function(Tree) {
  out <- Tree$tip.label[grep("PESP|AUAN", Tree$tip.label, perl=T)]
  res <- drop.tip(Tree, out)
  return(res)
}

# function to find the names of relevant tips 
find.tips <- function(tree, tips) {
  tips.pattern <- paste(tips, collapse="|")
  tips.in.tree <- tree$tip.label[grep(pattern = tips.pattern, x=tree$tip.label, perl = TRUE)]
  tips.in.tree
}

# function to extract clades of minimum size when lineage is not monophyletic
extractLargestCladesWhenNotMonophyletic <- function(phy, group, min_taxa=3, plot_trees=TRUE) {
  # all trees are returned as class 'multiPhylo' for easier combining afterwards
  if (is.monophyletic(phy = phy, tips = group)) {
    final <- drop.tip(phy, tip = phy$tip.label[!(phy$tip.label %in% group)])
    class(final) <- "phylo"
    if (plot_trees) {
      plot(final, cex=.6)
    }
  } else {
    # get all clades
    pp.phy <- prop.part(phy)
    # rename clades with their labels
    named.pp.phy <- lapply(1:length(pp.phy), function(x) {y <- attr(pp.phy, 'labels')[unlist(pp.phy[x])]; return(y)})
    # find clades whose tips are in the wanted tip vector
    completely.nested.in.wanted <- sapply(named.pp.phy, function(x) all(x %in% group))
    # subset the list to find those clades (sets of tips) that are completely nested
    nested.clades <- named.pp.phy[completely.nested.in.wanted]
    # find the largest clades that don't share any tips
    largest.nested <- sapply(1:length(nested.clades), function(x) all(unlist(nested.clades[x]) %in% unlist(nested.clades[-x])))
    # get the subtrees
    largest.subtrees <- lapply(nested.clades[!largest.nested], function(x) drop.tip(phy, tip = phy$tip.label[!(phy$tip.label %in% unlist(x))]))
    index <- sapply(largest.subtrees, Ntip) >= min_taxa
    if (sum(index)==0) {
      final <- NULL 
      #stop("No trees larger than min_taxa were found")
      cat("No trees larger than min_taxa were found\n")
    } else if (sum(index)==1){
      final <- largest.subtrees[index][[1]]
      class(final) <- "phylo"
    } else if (sum(index)>1) {
      final <- largest.subtrees[index]
      class(final) <- "multiPhylo"
    }
    if(plot_trees & sum(index)>1) {
      par(mfrow = c(1:length(final)))
      lapply(final, plot, cex=.6)
    }
  }
  return(final)
}

# function to extract a group when all members form a monophyletic group
extract.mono <- function(tree, tips) {
  tips.pattern <- paste(tips, collapse="|")
  tips.in.tree <- tree$tip.label[grep(pattern = tips.pattern, x=tree$tip.label, perl = TRUE)]
  if (length(tips.in.tree) >= 4) {
    mono <- try(is.monophyletic(tree, tips = tips.in.tree, reroot = TRUE))
  } else {
    mono <- FALSE
  }
  if (!is.null(mono) && mono) {
    res <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% tips.in.tree)])
    return(res)
  }
}


# Species tree
ST <- read.tree("original-data/absolute-time/absolute-time-rerooted-RT_38_72h.contree.tre")
ST <- root(ST, "BOPA", resolve.root = TRUE)
ST <- ladderize(ST)
ST$node.label <- NULL
#plot(ST)

# Species trees for smaller clades
thals <- c("DIPS","THOC","SKMA","THPS","THWE","POPS")
thal.ST <- drop.tip(ST, ST$tip.label[!(ST$tip.label %in% thals)])
penns <- c("CRAM","SESP","PHTR","GYFA","FRCY","NISP","CYCL","EUNA","FRSP","THAN","DITE","STSP","ASGL","TAPO","STUN")
penn.ST <- drop.tip(ST, ST$tip.label[!(ST$tip.label %in% penns)])

# load gene trees
tr <- read.tree("original-data/gene_trees.bootstrapped.min30intaxa.tre")
# rename tips
tr2 <- mclapply(tr, function(x) {x$tip.label <- renametips(x$tip.label)$X3; x}, mc.cores=detectCores()-1)

# A. 
# set: 3.1 K
# bootstrap filter: none
# clade filter: none

# drop outgroups (species tree doesn't have them)
tr3 <- mclapply(tr2, dropOutgroup, mc.cores = detectCores()-1)
class(tr3) <- "multiPhylo"
# write gene trees to file
write.tree(tr3, file="output/cladeAll-bs00/gene-trees.newick")
# write species tree to file in same folder
write.tree(ST, file="output/cladeAll-bs00/species-tree.newick")

write.tree(tr3, file="output/cladeAll-bs00-rearr40/gene-trees.newick")
write.tree(ST, file="output/cladeAll-bs00-rearr40/species-tree.newick")
write.tree(tr3, file="output/cladeAll-bs00-rearr50/gene-trees.newick")
write.tree(ST, file="output/cladeAll-bs00-rearr50/species-tree.newick")
write.tree(tr3, file="output/cladeAll-bs00-rearr60/gene-trees.newick")
write.tree(ST, file="output/cladeAll-bs00-rearr60/species-tree.newick")
write.tree(tr3, file="output/cladeAll-bs00-rearr70/gene-trees.newick")
write.tree(ST, file="output/cladeAll-bs00-rearr70/species-tree.newick")

# B. 
# set: 3.1 K
# bootstrap filter: 40
# clade filter: none

# filter such that mean BS >= 40
Index <- sapply(tr2, function(x) mean(as.numeric(x$node.label), na.rm=TRUE) >=40)
tr3 <- tr2[Index]
# drop outgroups (species tree doesn't have them)
tr4 <- mclapply(tr3, dropOutgroup, mc.cores = detectCores()-1)
class(tr4) <- "multiPhylo"
write.tree(tr4, file="output/cladeAll-bs40/gene-trees.newick")
write.tree(ST, file="output/cladeAll-bs40/species-tree.newick")

# C. 
# set: 3.1 K
# bootstrap filter: 50
# clade filter: none

# filter such that mean BS >= 50
Index <- sapply(tr2, function(x) mean(as.numeric(x$node.label), na.rm=TRUE) >=50)
tr3 <- tr2[Index]
# drop outgroups (species tree doesn't have them)
tr4 <- mclapply(tr3, dropOutgroup, mc.cores = detectCores()-1)
class(tr4) <- "multiPhylo"
write.tree(tr4, file="output/cladeAll-bs50/gene-trees.newick")
write.tree(ST, file="output/cladeAll-bs50/species-tree.newick")

#####
# extract clades of minimum size when lineage is not monophyletic

# thals
# find tips for each trees
tips.in.each.tree <- lapply(tr2, function(x) find.tips(tree = x, tips = thals))
length(tr2)==length(tips.in.each.tree)

clades.in.each.tree <- lapply(1:length(tr2), function(x)
  try(extractLargestCladesWhenNotMonophyletic(phy = tr2[[x]], group = tips.in.each.tree[[x]], min_taxa = 4, plot_trees = FALSE))   )
Index <- sapply(clades.in.each.tree, function(x) class(x) %in% c("phylo", "multiPhylo"))
sum(Index)
final.clades.in.each.tree <- clades.in.each.tree[Index]
flat.list.of.trees <- do.call(c, final.clades.in.each.tree)
length(flat.list.of.trees)

# K. 
# set 3.1 K
# bootstrap filter: 00
# clade filter: thalassiosirales

write.tree(flat.list.of.trees, file = "output/cladeThals-bs00/gene-trees.newick")
write.tree(thal.ST, file="output/cladeThals-bs00/species-tree.newick")

write.tree(flat.list.of.trees, file = "output/cladeThals-bs00-rearr40/gene-trees.newick")
write.tree(thal.ST, file="output/cladeThals-bs00-rearr40/species-tree.newick")
write.tree(flat.list.of.trees, file = "output/cladeThals-bs00-rearr50/gene-trees.newick")
write.tree(thal.ST, file="output/cladeThals-bs00-rearr50/species-tree.newick")
write.tree(flat.list.of.trees, file = "output/cladeThals-bs00-rearr60/gene-trees.newick")
write.tree(thal.ST, file="output/cladeThals-bs00-rearr60/species-tree.newick")
write.tree(flat.list.of.trees, file = "output/cladeThals-bs00-rearr70/gene-trees.newick")
write.tree(thal.ST, file="output/cladeThals-bs00-rearr70/species-tree.newick")

# L. 
# set 3.1 K
# bootstrap filter: 40
# clade filter: thalassiosirales

Index <- sapply(flat.list.of.trees, function(x) mean(as.numeric(x$node.label), na.rm=TRUE) >=40)
thal.mono.trees.2 <- flat.list.of.trees[Index]
write.tree(thal.mono.trees.2, file = "output/cladeThals-bs40/gene-trees.newick")
write.tree(thal.ST, file="output/cladeThals-bs40/species-tree.newick")

# M. 
# set 3.1 K
# bootstrap filter: 50
# clade filter: thalassiosirales

Index <- sapply(flat.list.of.trees, function(x) mean(as.numeric(x$node.label), na.rm=TRUE) >=50)
thal.mono.trees.2 <- flat.list.of.trees[Index]
write.tree(thal.mono.trees.2, file = "output/cladeThals-bs50/gene-trees.newick")
write.tree(thal.ST, file="output/cladeThals-bs50/species-tree.newick")

# N. 
# set 3.1 K
# bootstrap filter: 70
# clade filter: thalassiosirales

Index <- sapply(flat.list.of.trees, function(x) mean(as.numeric(x$node.label), na.rm=TRUE) >=70)
thal.mono.trees.2 <- flat.list.of.trees[Index]
write.tree(thal.mono.trees.2, file = "output/cladeThals-bs70/gene-trees.newick")
write.tree(thal.ST, file="output/cladeThals-bs70/species-tree.newick")

#####
# extract clades of minimum size when not monophyletic

# pennates
# find tips for each trees
tips.in.each.tree <- lapply(tr2, function(x) find.tips(tree = x, tips = penns))
length(tr2)==length(tips.in.each.tree)

clades.in.each.tree <- lapply(1:length(tr2), function(x)
  try(extractLargestCladesWhenNotMonophyletic(phy = tr2[[x]], group = tips.in.each.tree[[x]], min_taxa = 4, plot_trees = FALSE))   )
Index <- sapply(clades.in.each.tree, function(x) class(x) %in% c("phylo", "multiPhylo"))
sum(Index)
final.clades.in.each.tree <- clades.in.each.tree[Index]
flat.list.of.trees <- do.call(c, final.clades.in.each.tree)
length(flat.list.of.trees)

# O. 
# set 3.1 K
# bootstrap filter: 00
# clade filter: pennates

write.tree(flat.list.of.trees, file = "output/cladePenns-bs00/gene-trees.newick")
write.tree(penn.ST, file="output/cladePenns-bs00/species-tree.newick")

write.tree(flat.list.of.trees, file = "output/cladePenns-bs00-rearr40/gene-trees.newick")
write.tree(penn.ST, file="output/cladePenns-bs00-rearr40/species-tree.newick")
write.tree(flat.list.of.trees, file = "output/cladePenns-bs00-rearr50/gene-trees.newick")
write.tree(penn.ST, file="output/cladePenns-bs00-rearr50/species-tree.newick")
write.tree(flat.list.of.trees, file = "output/cladePenns-bs00-rearr60/gene-trees.newick")
write.tree(penn.ST, file="output/cladePenns-bs00-rearr60/species-tree.newick")
write.tree(flat.list.of.trees, file = "output/cladePenns-bs00-rearr70/gene-trees.newick")
write.tree(penn.ST, file="output/cladePenns-bs00-rearr70/species-tree.newick")


# P. 
# set 3.1 K
# bootstrap filter: 40
# clade filter: pennates

Index <- sapply(flat.list.of.trees, function(x) mean(as.numeric(x$node.label), na.rm=TRUE) >=40)
thal.mono.trees.2 <- flat.list.of.trees[Index]
write.tree(thal.mono.trees.2, file = "output/cladePenns-bs40/gene-trees.newick")
write.tree(penn.ST, file="output/cladePenns-bs40/species-tree.newick")

# Q. 
# set 3.1 K
# bootstrap filter: 50
# clade filter: pennates

Index <- sapply(flat.list.of.trees, function(x) mean(as.numeric(x$node.label), na.rm=TRUE) >=50)
thal.mono.trees.2 <- flat.list.of.trees[Index]
write.tree(thal.mono.trees.2, file = "output/cladePenns-bs50/gene-trees.newick")
write.tree(penn.ST, file="output/cladePenns-bs50/species-tree.newick")

# R. 
# set 3.1 K
# bootstrap filter: 70
# clade filter: pennates

Index <- sapply(flat.list.of.trees, function(x) mean(as.numeric(x$node.label), na.rm=TRUE) >=70)
thal.mono.trees.2 <- flat.list.of.trees[Index]
write.tree(thal.mono.trees.2, file = "output/cladePenns-bs70/gene-trees.newick")
write.tree(penn.ST, file="output/cladePenns-bs70/species-tree.newick")


#####
## NOT USED 

## EXTRACTING CLADES ONLY WHEN MONOPHYLETIC

# # extract thalassiosirales clades
# # only when tips are monophyletic
# # with this filter we look for WGD within the group
# 
# # find gene trees where thals are monophyletic
# thals <- c("DIPS","THOC","SKMA","THPS","THWE","POPS")
# 
# 
# thal.mono.trees <- lapply(tr3, extract.mono, thals)
# index <- sapply(thal.mono.trees, is.null)
# thal.mono.trees <- thal.mono.trees[!index]
# class(thal.mono.trees) <- "multiPhylo"
# thal.mono.trees
# 
# # D. 
# # set 3.1 K
# # bootstrap filter: none
# # clade filter: thalassiosirales
# 
# # just write above trees to file
# write.tree(thal.mono.trees, file="output/set3.1_bs00_cladeThals/gene-trees.newick")
# write.tree(thal.ST, file="output/set3.1_bs00_cladeThals/species-tree.newick")
# 
# # E. 
# # set 3.1 K
# # bootstrap filter: 40
# # clade filter: thalassiosirales
# 
# Index <- sapply(thal.mono.trees, function(x) mean(as.numeric(x$node.label), na.rm=TRUE) >=40)
# thal.mono.trees.2 <- thal.mono.trees[Index]
# write.tree(thal.mono.trees.2, file="output/set3.1_bs40_cladeThals/gene-trees.newick")
# write.tree(thal.ST, file="output/set3.1_bs40_cladeThals/species-tree.newick")
# 
# # F. 
# # set 3.1 K
# # bootstrap filter: 50
# # clade filter: thalassiosirales
# 
# Index <- sapply(thal.mono.trees, function(x) mean(as.numeric(x$node.label), na.rm=TRUE) >=50)
# thal.mono.trees.2 <- thal.mono.trees[Index]
# write.tree(thal.mono.trees.2, file="output/set3.1_bs50_cladeThals/gene-trees.newick")
# write.tree(thal.ST, file="output/set3.1_bs50_cladeThals/species-tree.newick")
# 
# # G. 
# # set 3.1 K
# # bootstrap filter: 70
# # clade filter: thalassiosirales
# 
# Index <- sapply(thal.mono.trees, function(x) mean(as.numeric(x$node.label), na.rm=TRUE) >=70)
# thal.mono.trees.2 <- thal.mono.trees[Index]
# write.tree(thal.mono.trees.2, file="output/set3.1_bs70_cladeThals/gene-trees.newick")
# write.tree(thal.ST, file="output/set3.1_bs70_cladeThals/species-tree.newick")
# 
# # extract pennate clades
# # only when tips are monophyletic
# # with this filter we look for WGD within the group
# 
# # find gene trees where pennates are monophyletic
# penns <- c("CRAM","SESP","PHTR","GYFA","FRCY","NISP","CYCL","EUNA","FRSP","THAN","DITE","STSP","ASGL","TAPO","STUN")
# penn.ST <- drop.tip(ST, ST$tip.label[!(ST$tip.label %in% penns)])
# 
# penns.mono.trees <- lapply(tr3, extract.mono, penns)
# index <- sapply(penns.mono.trees, is.null)
# penns.mono.trees <- penns.mono.trees[!index]
# class(penns.mono.trees) <- "multiPhylo"
# penns.mono.trees
# 
# # H. 
# # set 3.1 K
# # bootstrap filter: none
# # clade filter: pennates (sin Attheya)
# 
# # just write above trees to file
# write.tree(penns.mono.trees, file="output/set3.1_bs00_cladePenns/gene-trees.newick")
# write.tree(penn.ST, file="output/set3.1_bs00_cladePenns/species-tree.newick")
# 
# # I. 
# # set 3.1 K
# # bootstrap filter: 40
# # clade filter: pennates
# 
# Index <- sapply(penns.mono.trees, function(x) mean(as.numeric(x$node.label), na.rm=TRUE) >=40)
# penns.mono.trees.2 <- penns.mono.trees[Index]
# write.tree(penns.mono.trees.2, file="output/set3.1_bs40_cladePenns/gene-trees.newick")
# write.tree(penn.ST, file="output/set3.1_bs40_cladePenns/species-tree.newick")
# 
# # J. 
# # set 3.1 K
# # bootstrap filter: 50
# # clade filter: pennates
# 
# Index <- sapply(penns.mono.trees, function(x) mean(as.numeric(x$node.label), na.rm=TRUE) >=50)
# penns.mono.trees.2 <- penns.mono.trees[Index]
# write.tree(penns.mono.trees.2, file="output/set3.1_bs50_cladePenns/gene-trees.newick")
# write.tree(penn.ST, file="output/set3.1_bs50_cladePenns/species-tree.newick")
# 
