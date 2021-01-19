# joy plot for Ks data

library("ggjoy")
library("ape")
library("ggtree")
library("tidyverse")
library("viridis")

rm(list=ls())

ks.joy.plot <- function(PathToFiles, Title, Pipeline=c("Johnson", "McKain"), Legend=FALSE) {
  Files <- list.files(path = PathToFiles, pattern = "txt", all.files = TRUE, full.names = TRUE, recursive = TRUE)
  if (Pipeline=="Johnson") {
    Names <- strsplit(Files, "/") %>% do.call(rbind, .) %>% .[,5] %>% gsub(".kaks.txt", "", .)
  } else {
    Names <- strsplit(Files, "/") %>% do.call(rbind, .) %>% .[,5] %>% strsplit(., "_") %>% do.call(rbind, .) %>% .[,1]
  }
  print(cat(Names))
  TaxGroups <- read.csv("original-data/taxa-groups.csv", header=F)
  colnames(TaxGroups) <- c("Species", "Group")
  
  species.tree <- read.tree("original-data/absolute-time/absolute-time-rerooted-RT_38_72h.contree.tre")
  TaxGroups <- TaxGroups[match(species.tree$tip.label, TaxGroups$Species),]
   
  Ks.list <- list()
  for (i in seq_along(Files)) {
    Ks.list[[i]] <- read.table(Files[i], header=TRUE)
    Ks.list[[i]]$Species <- Names[i]
  }
  
  Ks.data <- Ks.list %>% bind_rows() %>% filter(!(Species =="AUAN")) %>% 
    inner_join(TaxGroups, ., by="Species") %>% 
    mutate(Group2=factor(Group, levels=c("Raphid", "Araphid", "Thalassiosirales", "Polar", "Radial"))) %>%
    mutate(Species2=factor(Species, levels=species.tree$tip.label))
  print(head(Ks.data))
  
  P <- ggplot(Ks.data, aes(x=ks_value, y=Species2, fill=Group2)) + 
    geom_joy(scale=4, size=.5, colour="grey23", alpha=.2, panel_scaling = FALSE) +
    theme_joy(grid = TRUE, font_size = 8, line_size = .5) +
    scale_fill_brewer(palette = "Set1", type = "qual", name="") +
    scale_x_continuous(expand = c(0.01, 0)) +
    scale_y_discrete(expand = c(0.01, 0)) +
    xlab("Ds") + ylab("") + #xlim(0,2) +
    ggtitle(label = "", subtitle = Title)
  if(Legend) {P <- P+theme(legend.position="top")} else {P <- P+theme(legend.position="n")}
  # +
  #   facet_wrap("Group")
  return(P)
  }

jgy <- ks.joy.plot(PathToFiles = "original-data/Ks-data/Ks_values_full/johnson.gy", Title = "Johnson GY", Pipeline = "Johnson")
jgn <- ks.joy.plot(PathToFiles = "original-data/Ks-data/Ks_values_full/johnson.yn", Title = "Johnson YN", Pipeline = "Johnson")
mba <- ks.joy.plot(PathToFiles = "original-data/Ks-data/Ks_values_full/mckain_blast_all", Title = "McKain blast all", Pipeline = "McKain")
mgc <- ks.joy.plot(PathToFiles = "original-data/Ks-data/Ks_values_full/mckain_blast_gene_collapsed", Title = "McKain blast gene collapsed", Pipeline = "McKain")

#ggsave(gridExtra::grid.arrange(jgy, jgn, mba, mgc, ncol=4), filename="Ks.joy.plot.pdf", path="figures", width=12, height=8)


### combine for same scale

loadData <- function(PathToFiles, Title, Pipeline=c("Johnson", "McKain"), ReturnList=FALSE) {
  Files <- list.files(path = PathToFiles, pattern = "txt", all.files = TRUE, full.names = TRUE, recursive = TRUE)
  if (Pipeline=="Johnson") {
    Names <- strsplit(Files, "/") %>% do.call(rbind, .) %>% .[,5] %>% gsub(".kaks.txt", "", .)
  } else {
    Names <- strsplit(Files, "/") %>% do.call(rbind, .) %>% .[,5] %>% strsplit(., "_") %>% do.call(rbind, .) %>% .[,1]
  }
  print(cat(Names))
  TaxGroups <- read.csv("original-data/taxa-groups.csv", header=F)
  colnames(TaxGroups) <- c("Species", "Group")
  
  species.tree <- read.tree("original-data/absolute-time/absolute-time-rerooted-RT_38_72h.contree.tre")
  species.tree <- ladderize(species.tree)
  TaxGroups <- TaxGroups[match(species.tree$tip.label, TaxGroups$Species),]
  
  Ks.list <- list()
  for (i in seq_along(Files)) {
    Ks.list[[i]] <- read.table(Files[i], header=TRUE)
    Ks.list[[i]]$Species <- Names[i]
  }
  Ks.data <- Ks.list %>% bind_rows() %>% filter(!(Species =="AUAN")) %>% 
    inner_join(TaxGroups, ., by="Species") %>% 
    mutate(Group2=factor(Group, levels=c("Raphid", "Araphid", "Thalassiosirales", "Polar", "Radial"))) %>%
    mutate(Species2=factor(Species, levels=species.tree$tip.label)) %>%
    mutate(Method=Title)
  
  if (ReturnList) {
    return(list(Ks.data, Ks.list))
  } else {
    return(Ks.data)
  }
}

A <- loadData(PathToFiles = "original-data/Ks-data/Ks_values_full/johnson.gy", Title = "Johnson GY", Pipeline = "Johnson")
B <- loadData(PathToFiles = "original-data/Ks-data/Ks_values_full/johnson.yn", Title = "Johnson YN", Pipeline = "Johnson")
C <- loadData(PathToFiles = "original-data/Ks-data/Ks_values_full/mckain_blast_all", Title = "McKain blast all", Pipeline = "McKain")
D <- loadData(PathToFiles = "original-data/Ks-data/Ks_values_full/mckain_blast_gene_collapsed", Title = "McKain blast gene collapsed", Pipeline = "McKain")

EE <- bind_rows(A,B,C,D)

P <- ggplot(EE, aes(x=ks_value, y=Species2, fill=Group2)) + 
  #geom_joy(stat = "binline", bins = 50, scale = 0.95, draw_baseline = FALSE) +
  geom_joy(scale=4, size=.3, colour="grey23", alpha=.4, panel_scaling = FALSE) +
  theme_joy(grid = TRUE, font_size = 8, line_size = .5) +
  scale_fill_brewer(palette = "Set1", type = "qual", name="") +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(expand = c(0.01, 0)) +
  xlab("Ds") + ylab("") + #xlim(NA,2) +
  ggtitle(label = "") +
  facet_wrap(~Method, ncol=4) +
  theme(legend.position = "top")
P
#ggsave(P, filename="Ks.joy.plot.same.scale.pdf", path="figures", width=8, height=8)
A$Species3 <- factor(A$Species2, levels=c("BOPA", "COHY", "LEDA", "PRAL", "AUSU", "STTU",
                                          "ACSU", "RHSE", "LASH", "MIPO", "ODAU", "EUAN", 
                                          "HESI", "CHAF", "CHSP", "DIBR", "POPS", "THWE",
                                          "THPS", "SKMA", "DIPS", "THOC", "ATSE", "STUN",
                                          "ASGL", "TAPO", "STSP", "DITE", "FRSP", "THAN",
                                          "EUNA", "CYCL", "FRCY", "NISP", "GYFA", "PHTR", "CRAM", "SESP"))


P <- ggplot(A, aes(x=ks_value, y=Species3, fill=Species3)) + 
  #geom_joy(stat = "binline", bins = 50, scale = 0.95, draw_baseline = FALSE) +
  geom_joy(scale=4, size=.3, colour="grey23", alpha=.75, panel_scaling = FALSE) +
  theme_classic() +
  #theme_joy(grid = TRUE, font_size = 8, line_size = .5) +
  #scale_fill_brewer(palette = "Set1", type = "qual", name="") +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(expand = c(0.01, 0)) +
  xlab("Ks") + ylab("") + #xlim(NA,2) +
  ggtitle(label = "") +
  #facet_wrap(~Method, ncol=4) +
  theme(legend.position = "none") +
  #scale_fill_cyclical(values = c("black", "white")) 
  scale_fill_cyclical(values = c("#756bb1", "#bcbddc")) 
P

ggsave(P, filename="Ks.john.gy.reorder.pdf", path="figures", width=6, height=11.1)

# just Thals

JT <- EE %>% filter(Method=="Johnson GY") %>% filter(Group2=="Thalassiosirales")

P2 <- ggplot(JT, aes(x=ks_value, y=Species2, fill=Group2, height=..density..)) + 
  geom_joy(stat = "binline", bins = 10, size=.3, alpha=.4, scale = 0.95, draw_baseline = FALSE) +
  #geom_joy(scale=3, size=.3, colour="grey23", alpha=.4, panel_scaling = FALSE, fill="grey60") +
  theme_joy(grid = TRUE, font_size = 8, line_size = .5) +
  #scale_fill_brewer(palette = "Set1", type = "qual", name="") +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(expand = c(0.01, 0)) +
  xlab("Ks") + ylab("") + #xlim(NA,2) +
  ggtitle(label = "") +
#   facet_wrap(~Method, ncol=4) +
  theme(legend.position = "none")
P2

ggsave(P2, filename="Thals.joy.pdf", path="figures", width=8, height=8)


# faceted standard histograms/densities
# with peaks

glimpse(A)
#install.packages("mclust")
library("mclust")
A <- loadData(PathToFiles = "Ks-data/Ks_values_full/johnson.gy", Title = "Johnson GY", Pipeline = "Johnson", ReturnList = TRUE)
A2 <- A[[2]]
XX <- lapply(A2, function(x) densityMclust(x$ks_value))
plot(XX[[1]], what = "density", data = A2[[1]]$ks_value, breaks = 15)

Dens <- ggplot(A, aes(x=ks_value, group=Species2)) +
  geom_density2d(stat = "density", colour="grey23") +
  facet_wrap("Species2", strip.position = "top", ncol=8) +
#  scale_colour_brewer(palette = "Set1", type = "qual", name="") +
  ylim(0,3) +
  theme_minimal() +
  theme(legend.position="top") 
Dens 
ggsave(Dens, filename = "Ks.density.pdf", path = "figures", width = 10, height=10)

Hist <- ggplot(A, aes(x=ks_value, group=Species2)) +
  geom_histogram( colour="white") +
  facet_wrap("Species2", strip.position = "top", ncol=8) +
  #  scale_colour_brewer(palette = "Set1", type = "qual", name="") +
  ylim(0,500) +
  theme_minimal() +
  theme(legend.position="top") 
Hist

#### some made up data on peaks

P2 <- ggplot(JT, aes(x=ks_value, y=Species2, fill=Species2, height=..density..)) + 
  geom_joy(stat = "binline", bins = 30, size=.5, alpha=.4, scale = 1, draw_baseline = FALSE) +
  #geom_joy(scale=3, size=.3, colour="grey23", alpha=.4, panel_scaling = FALSE, fill="grey60") +
  #theme_joy(grid = TRUE, font_size = 8, line_size = .5) +
  #scale_fill_brewer(palette = "Set1", type = "qual", name="") +
  facet_wrap(~Species, ncol=1, scales="free_y", as.table = FALSE, strip.position = "left") +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(expand = c(0.01, 0)) +
  xlab("Ks") + ylab("") + #xlim(NA,2) +
  ggtitle(label = "") +
  theme_minimal() +
  theme(legend.position = "none", 
        #axis.text.y = element_text(angle=90, hjust = -.75)
        axis.text.y = element_blank()) +
  scale_fill_cyclical(values = c("grey20", "grey90")) 
P2

# P3 <- ggplot(JT, aes(x=ks_value, group=Species2, height=..density..)) +
#   #geom_histogram(bins=50, colour="white") +
#   geom_density(colour="white", fill="gray23") +
#   facet_wrap(~Species, ncol=1, as.table = FALSE, strip.position = "left") +
#   # scale_x_continuous(expand = c(0.01, 0)) +
#   # scale_y_discrete(expand = c(0.01, 0)) +
#   xlab("Ks") + ylab("") + #xlim(NA,2) +
#   ggtitle(label = "") +
#   theme_minimal() +
#   theme(legend.position = "none") 
#         #axis.text.y = element_text(angle=90, hjust = -.75)
#         #axis.text.y = element_blank())
# P3

# ## Bubbles (fake data)
# spec <- unique(JT$Species)
# peak <- rep(1:3, length(spec))
# Pipeline <- c("JGY", "JYN", "MBA", "MGC")
# Species <- expand.grid(spec, Pipeline) %>% mutate(Comb=paste(Var1, Var2, sep=".")) %>% arrange(Var1)
# peakInfo2 <- data.frame(Species=rep(Species$Var1, 3),
#                         SpecPipe=rep(Species$Comb, 3),
#                         Peak=rep(peak,each=4), 
#                         Ks=runif(min = 0.5, max=2, n=nrow(Species)*3), 
#                         Genes=sample(200:1000, nrow(Species)*3), 
#                         Pipeline=rep(Pipeline, length(peak)),
#                         stringsAsFactors = FALSE)
# peakInfo2$Species <- factor(peakInfo2$Species, levels=c("DIPS", "THOC", "SKMA", "THPS", "THWE", "POPS", "DIBR"))
# 
# # remove some rows to mimick missing peaks
# peakInfo2[sample(1:nrow(peakInfo2), 40),4:5] <- NA

## the tree
tr <- read.tree("original-data/absolute-time/absolute-time-rerooted-RT_38_72h.contree.tre")
tr <- ladderize(tr)
# plot
plot.tr <- ggtree(tr=tr) + 
  geom_tiplab() + 
  theme_tree2()
plot.tr

## Bubbles

Kstab <- read.csv("output/Ks_table.txt") %>% 
  filter(!Species=="AUAN") %>% 
  # filter(Total_paralog_pairs>) %>%
  filter(Paralog_pairs_under_peak/Total_paralog_pairs>0.2) %>% 
  mutate(SpecPipe=paste(Species, Pipeline, sep="|"))
#Kstab <- Kstab %>% mutate(Species2=factor(Species, levels=rev(tr$tip.label)))
Kstab <- Kstab %>% mutate(Species2=factor(Species, levels=c("BOPA", "COHY", "LEDA", "PRAL", "AUSU", "STTU",
                                                            "ACSU", "RHSE", "LASH", "MIPO", "ODAU", "EUAN", 
                                                            "HESI", "CHAF", "CHSP", "DIBR", "POPS", "THWE",
                                                            "THPS", "SKMA", "DIPS", "THOC", "ATSE", "STUN",
                                                            "ASGL", "TAPO", "STSP", "DITE", "FRSP", "THAN",
                                                            "EUNA", "CYCL", "FRCY", "NISP", "GYFA", "PHTR", "CRAM", "SESP")))

Bubbles <- ggplot(data=Kstab, aes(group=Species2, y=Species2, x=Ks_value, colour=Pipeline, size=Paralog_pairs_under_peak)) +
  geom_point(alpha=.8, stroke=1, pch=24) + 
  #facet_wrap(~Species, ncol=1, scales="free_y", as.table = FALSE, strip.position = "left") +
  theme_minimal() +
  xlab("Ks") +
  ylab("") +
  scale_colour_viridis(discrete=TRUE) +
  scale_size_continuous(name="Paralog pairs\nunder peak", breaks = c(300, 600, 1000), labels = c(300, 600, 1000), range = c(0.5, 7)) +
  scale_x_continuous(breaks = seq(0,2,.5), labels = seq(0,2,.5), limits = c(0,2)) +
  theme(#legend.position="top",
        legend.box = "vertical",
        legend.box.just = "left",
        legend.spacing.y = unit(.05, "cm")
    # axis.text.y = element_blank()
        )
Bubbles

ggsave(Bubbles, filename = "paralogs-under-peak.pdf", path = "figures", width = 5, height=12)
YY <- cowplot::plot_grid(plot.tr, Bubbles, labels = c("A", "B"), align = "hv", axis = "b", rel_widths = c(1,2))
ggsave(YY, filename = "Tree-Ks-peaks.pdf", path = "figures", width = 9, height=7)

#### bubbles facet

Bubbles2 <- ggplot(data=Kstab, aes(group=Species2, y=Pipeline, x=Ks_value, fill=Pipeline, size=Paralog_pairs_under_peak/Total_paralog_pairs*100)) +
  geom_point(alpha=.6, stroke=1, pch=21) + 
  facet_wrap(~Species, ncol=5, scales="free_y", as.table = TRUE, strip.position = "top") +
  theme_minimal() +
  xlab("Ks") +
  ylab("") +
  scale_fill_viridis(discrete=TRUE) +
  scale_size_continuous(name="Paralog pairs\nunder peak\n(% of total)", breaks = c(25, 45, 65), labels = c("25%", "45%", "65%"), range = c(2.5, 6.8)) +
  scale_x_continuous(breaks = seq(0,2,.5), labels = seq(0,2,.5), limits = c(0,2)) +
  theme(legend.position=c(0.8, 0.05),
    legend.box = "vertical",
    legend.box.just = "left",
    legend.direction = "horizontal",
    legend.spacing.y = unit(.05, "cm"),
    legend.key.size = unit(.5, "cm"),
    #axis.text.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=.5)
  ) +
  guides(fill = guide_legend(override.aes = list(size=4)))
#Bubbles2

ggsave(Bubbles2, filename = "Ks-peaks-faceted.pdf", path = "figures", width = 12, height=12)

### Tiles (not used)
spec <- unique(JT$Species)
peak <- rep(1:3, length(spec))
peaksJGY <- rep(c(0.5, 1,2.5), length(spec))
peaksJYN <- peaksJGY +.2
genesJGY <- sample(1:1000, length(peaksJGY))
genesJYN <- sample(200:1500, length(peaksJGY))

peakInfo <- data.frame(Species=rep(spec, each=3), peak,
                       peaksJGY, genesJGY, JGYtext=paste(peaksJGY, "\n(", genesJGY, ")", sep=""),
                       peaksJYN, genesJYN, JYNtext=paste(peaksJYN, "\n(", genesJYN, ")", sep=""),
                       stringsAsFactors = FALSE)
peakInfo$Species <- factor(peakInfo$Species, levels=c("DIPS", "THOC", "SKMA", "THPS", "THWE", "POPS", "DIBR"))
# remove some rows to mimick missing peaks
peakInfo[c(1,10,11,19:21),3:8] <- NA

peakInfoLong <- peakInfo %>% group_by(Species, peak) %>%
  select(-ends_with("text")) %>%
  gather(Case, Value, peaksJGY:genesJYN) %>%
  mutate(Case2=factor(Case, levels=c("peaksJGY", "genesJGY", "peaksJYN", "genesJYN") )) %>%
  mutate(Type=ifelse(Case2 %in% c("peaksJGY", "peaksJYN"), "peaks", "genes")) %>%
  mutate(Pipeline=ifelse(Case2 %in% c("genesJGY", "peaksJGY"), "JGY", "JYN")) %>%
  mutate(Value2=ifelse(Value> 2, Value/1000, Value))

peakInfoLong2 <- peakInfo %>% group_by(Species, peak) %>%
  select(-starts_with("peaks")) %>%
  select(-starts_with("genes")) %>%
  gather(Case, Value, JGYtext:JYNtext) %>%
  mutate(Pipeline=ifelse(Case=="JGYtext", "JGY", "JYN"))

Tiles <- ggplot(data=peakInfoLong, aes(group=Species, y=Species, x=peak)) +
  geom_tile(colour="white", size=1, aes(fill=Value), data=peakInfoLong, alpha=.8) +
  #geom_tile(colour="grey23", size=.3, fill="white", data=peakInfoLong) +
  geom_text(aes(group=Species, y=Species, x=peak, label=Value), size=3, colour="black",fontface="bold", data=peakInfoLong2) +
  facet_wrap(~Pipeline, ncol=2) +
  theme_minimal() +
  xlab("Ks Peak") +
  ylab("") +
  scale_fill_continuous(guide=FALSE, low = "white", high = "firebrick4")
Tiles
# 
# XX <- cowplot::plot_grid(P2, Tiles, align = "hv", axis = "b", ncol = 2, rel_widths = c(1,1.5))
# XX
# ggsave(XX, filename = "Ks-with-info.pdf", path = "figures", width = 7, height=5)
