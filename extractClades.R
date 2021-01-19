
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
