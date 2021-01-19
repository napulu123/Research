# amoebas Ks with ridges
setwd("~/Desktop/Research/ks/")

##directory contains ks files and species tree file in data directory
##ks_files extension ".kaks.txt"
##species_tree ".tre"

rm(list=ls())
library(tidyverse)
library(ggridges)
library(ape)
library(ggtree)
#installed.packages("cowplots")
library(cowplot)
#installed.packages("patchwork")
library(patchwork)

# list all ks files in your data directory
all_ks_files <- list.files(path = "data/", pattern = ".kaks.txt", full.names = TRUE)
# load them all into R
all_ks_tables <- map(all_ks_files, function(x) readr::read_tsv(x))
# make names
all_ks_names <- stringr::str_extract(all_ks_files, regex("^.*_")) %>% stringr::str_replace(., regex("data//"), "") %>% stringr::str_replace("_", "")
all_ks_names <- tolower(all_ks_names)
# assign names
names(all_ks_tables) <- all_ks_names
# bind the tables into a single data frame
all_ks_df <- dplyr::bind_rows(all_ks_tables, .id="Species")
# extract just the ka/ks ratio column, if you wan't to
ks_simple <- all_ks_df %>% dplyr::select(Species, Ks)


##test code for adding clade information
##I got tired of trying so I stopped


# add a made up column for group or clade
# first way, might be cumbersome and error prone
# ks_simple$Clade <- case_when(
#   ks_simple$Species=="Mayocant" ~ "X",
#   ks_simple$Species=="Rhizlibe" ~ "X",
#   ks_simple$Species=="Rhizsaxo" ~ "Y",
#   ks_simple$Species=="RipDP13"  ~ "Z"
# )

# second way
# use a look up table
# make it in excel (or osmething) and load it here.
# then join
# species_lookup <- tibble(Species=all_ks_names, Clade=c("X", "X", "Y", "Z"))
# # join this table to the ks_simple table by Species
# ks_simple_w_clades <- dplyr::left_join(ks_simple, species_lookup, by="Species")



# add the tree
#load
tree <- read.tree("data/Amoebozoa_reroot.tre")
tree$tip.label <- tolower(tree$tip.label)
# rename some species on the tree
tree$tip.label[tree$tip.label == "acangen"] <- "acancast"
# tree$tip.label[tree$tip.label == "rhizlibe"] <- "acancast"
# tree$tip.label[tree$tip.label == "tricyna"] <- "tricyana"

# check if the names ks tables are all in the names of tips in the tree
sp_not_on_tree <- all_ks_names[!all_ks_names %in% tree$tip.label]
# remove those two species from the ks data
ks_simple_reduced <- ks_simple %>% filter(!Species %in% sp_not_on_tree)
new_ks_names <- ks_simple_reduced %>% select(Species) %>% distinct() %>% unlist %>% unname %>% tolower

# check again
new_ks_names[!new_ks_names %in% tree$tip.label]

# drop tips from the tree if they are not in the Ks data
# find the unwanted tips 
tree_tips_not_in_ks_data <- tree$tip.label[!tree$tip.label %in% new_ks_names]
# drop them from the tree
tree_reduced <- drop.tip(tree, tree_tips_not_in_ks_data)
# smooth branches just for plotting
tree_smooth <- chronopl(tree_reduced, lambda=0.001)

# Tree with Ks ridges

# Make the tree
tree_plot <- ggtree(tree_smooth) #+ theme_tree2() #+ geom_tiplab(size=3)
tree_levels <- tree_plot$data %>% dplyr::select(label, y) %>% drop_na() %>% arrange(y) %>% select(label) %>% unlist %>% unname

# make ggridges plot, same as above except reorder by the tree coordinates
ks_reordered <- ks_simple_reduced
ks_reordered$Species <- factor(ks_reordered$Species, level=tree_levels)

##change   scale_x_continuous(limits=c(NA,100)) 
## to   scale_x_continuous(limits=c(50,100)) to see speaks zoomed in

##make ggridges plot
ridges <- ggplot(data=ks_reordered, aes(x=Ks, y=Species)) +
  geom_density_ridges(scale=5, size=.5, colour="grey23", alpha=.5, panel_scaling = FALSE) +
  scale_x_continuous(limits=c(NA,100)) +
  theme(axis.text.y = element_text(size=8)) +
  labs(y="")
##generate Ks curve plots
ridges
# cowplot::plot_grid(tree_plot, ridges, ncol=2)
tree_plot + ridges


##adjust plot with inkscape, seems a bit off

