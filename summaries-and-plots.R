# compare notung with filters and rearrangements
# figure 3 in manuscript

library("tidyverse")
library("ape")

rm(list=ls())
dup.files <- list.files(path = "output", pattern = "duplication.txt", recursive = TRUE, full.names = TRUE)
dup.files <- dup.files[grep(pattern = "Thals|Penns", x = dup.files, perl = TRUE, invert = TRUE)]
dup.files <- dup.files[grep(pattern = "obsolete", x = dup.files, perl = TRUE, invert = TRUE)]
nnames <- strsplit(dup.files, "/") %>% unlist(., recursive=FALSE) %>% matrix(., nrow=7, byrow=TRUE) %>% .[,2] %>% gsub("cladeAll", "", .) %>% gsub("-", "", .)
dup.tables <- lapply(dup.files, function(x) read.table(x, header=T, row.names = 1, na.strings = "NaN", sep = "\t"))
for (i in 1:7) dup.tables[[i]]$Set <- nnames[i]
dup.tab.df <- bind_rows(dup.tables)

wanted <- data.frame(Names=paste("n", c(68,62,58,52,26,22), sep=""), 
                     Nodes= c("Diatoms\n- Corethron\n- Leptocylindrus", 
                              "Pennates\n+ Mediophytes\n+ Coscinodiscoids",
                              "Pennates\n+ Mediophytes", 
                              "Thalassiosirales\n- Porosira",
                              "Pennates\n- Striatella", 
                              "Pennates\n- Striatella\n- Asterionellopsis\n- Talaroneis"),
                     NodeNumbers= paste("Branch", LETTERS[1:6], sep=" "),
                     stringsAsFactors = FALSE)

dup.tab.df2 <- dup.tab.df %>% 
  group_by(Set) %>% 
  select(wanted$Names) %>%
  mutate_all(function(x) x > 0) %>% 
  mutate(N=n()) %>% 
  ungroup %>% 
  group_by(Set, N) %>% 
  summarise_at(.vars = vars(n68:n22), function(x) sum(x)) %>% 
  gather(Names, Count, n68:n22) %>% 
  mutate(Percent=Count/N*100)

# dup.tab.df3 <- 
#   dup.tab.df %>% 
#   group_by(Set) %>% 
#   #select(wanted$Names) %>%
#   mutate_all(function(x) x > 0) %>% 
#   mutate(N=n()) %>% 
#   ungroup %>% 
#   group_by(Set, N) %>% 
#   summarise_all(function(x) sum(x)) %>% 
#   gather(Names, Count, SESP:BOPA) %>% 
#   mutate(Percent=Count/N*100)


for.plot <- left_join(dup.tab.df2,wanted, by="Names") %>% 
  mutate(BS=case_when(Set == "bs00" ~ 0,
                      Set == "bs00rearr40" ~ 40,
                      Set == "bs00rearr50" ~ 50,
                      Set == "bs00rearr60" ~ 60,
                      Set == "bs00rearr70" ~ 70,
                      Set == "bs40" ~ 40,
                      Set == "bs50" ~ 50)) %>% 
  mutate(Type=ifelse(Set %in% c("bs00", "bs40", "bs50"), "Notung\nFiltered", "Notung\nRearranged"))

#### add the results from Yang pipeline
# 
# YangRes <- read.csv("output/Yang-duplication-counts-reordered-for-species-tree.csv", stringsAsFactors = FALSE, header=TRUE)
# YangRes <- filter(YangRes, sp.t.nodes %in% wanted$Name)
# YangRes2 <- data.frame(Set="Yang", N=NA, Names=wanted$Names, Count=NA, 
#                       Percent=YangRes$yang.dups,
#                       Nodes=wanted$Nodes, NodeNumbers=wanted$NodeNumbers, BS=40, 
#                       Type="Yang\nBS>40%", stringsAsFactors = FALSE)
# for.plot <- bind_rows(for.plot, YangRes2)

aa <- read.tree("output/RT_38_72h.contree.rooted.dup_mapped.fix.tre.dup_count.tre")
aa$node.label[aa$node.label==""] <- "0/1"
NL <- aa$node.label %>% strsplit(., "/") %>% do.call(rbind, .) %>% 
  data.frame(., stringsAsFactors=FALSE) %>% 
  mutate(Dups=as.numeric(X1), Tots=as.numeric(X2), Percent=Dups/Tots*100)
NL
# find the nodes on the tree
plot(aa, no.margin=T)
nodelabels(pch=21, cex=NL$Percent/10)
nodelabels(text=round(NL$Percent,2), font=2, cex=.7, bg="white", frame="n", adj=c(1.5, 1))
nodelabels(text=paste(NL$Dups, NL$Tots, sep="/"), font=2, cex=.5, bg="white", frame="n", adj=c(1.5, 1))
# extract the nodes
yang.wanted <- NL %>% filter(round(Percent,2) %in% c(43.8, 13.98, 53.18, 10.08, 25.63, 10.13))
yang.wanted <- yang.wanted %>% 
  mutate(Names=paste("n", c(68,62,58,52,26,22), sep=""), 
         Nodes= c("Diatoms\n- Corethron\n- Leptocylindrus", 
                  "Pennates\n+ Mediophytes\n+ Coscinodiscoids",
                  "Pennates\n+ Mediophytes", 
                  "Thalassiosirales\n- Porosira",
                  "Pennates\n- Striatella", 
                  "Pennates\n- Striatella\n- Asterionellopsis\n- Talaroneis"),
         NodeNumbers= paste("Branch", LETTERS[1:6], sep=" ")) %>% 
  select(-X1,X2) %>% 
  rename(N=Tots, Count=Dups) %>% 
  mutate(Set="Yang", Type="Yang\nBS>40%", BS=40) %>% 
  select(Set, N, Names, Count, Percent, Nodes, NodeNumbers, BS, Type)
yang.wanted

for.plot <- bind_rows(for.plot, yang.wanted)

# duplicate the default set with Type=Rearranged
def.set <- filter(for.plot, Set=="bs00")
def.set$Type <- "Notung\nRearranged"

for.plot <- bind_rows(for.plot, def.set)
for.plot$Type[for.plot$Type %in% c("Notung\nFiltered", "Notung\nRearranged") & for.plot$Set=="bs00"] <- "Notung\nOriginal"

# add a column for number of trees
# for.plot <- for.plot %>%
#   mutate(Ntrees=case_when(Set=="bs00" ~ 3163,
#                           Set=="bs00rearr40" ~ 3163,
#                           Set=="bs00rearr50" ~ 3163,
#                           Set=="bs00rearr60" ~ 3163,
#                           Set=="bs00rearr70" ~ 3163,
#                           Set=="bs40" ~ 1054,
#                           Set=="bs50" ~ 301,
#                           Set=="Yang" ~ 706))

for.plot$Type2 <- factor(for.plot$Type, levels=c("Notung\nOriginal", "Notung\nFiltered", "Notung\nRearranged", "Yang\nBS>40%"))
for.plot2 <- for.plot
for.plot2$Type[for.plot2$Type %in% c("Notung\nOriginal") & for.plot2$Set=="bs00"][1:6] <- "Notung\nFiltered"
for.plot2$Type[for.plot2$Type %in% c("Notung\nOriginal") & for.plot2$Set=="bs00"] <- "Notung\nRearranged"

XX <- ggplot(for.plot, aes(x=BS, y=Percent, group=Type2, colour=Type2, size=N)) + 
  geom_line(colour='grey23', size=.2, data=for.plot2, aes(x=BS, y=Percent, group=Type, linetype=Type)) + 
  geom_point(pch=21, alpha=.7, stroke=1, position = position_jitter(width=0)) + 
  scale_linetype(name="", guide=FALSE) +
  scale_colour_brewer(name="", type = "div", palette = "Set1") +
  scale_size(name="", range = c(0.5,4), breaks=c(500, 1000, 3000), labels=c("500\ntrees", "2000\ntrees", "3000\ntrees")) +
  #scale_fill_brewer(name="", type = "div", palette = "Set1") +
  facet_wrap("NodeNumbers", ncol=3, scales="free_x") +  
  theme_minimal() + 
  theme(legend.position = "top") +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  scale_y_continuous(name="Gene families with duplication (%)", breaks = seq(0,100,10), minor_breaks = seq(0,100,10), labels = seq(0,100,10), limits = c(-10,100), sec.axis = dup_axis()) +
  scale_x_continuous(name="Bootstrap threshold (%)", breaks = seq(0,100,20), minor_breaks = seq(0,100,10), limits = c(-10,80)) +
   labs(title="Figure 3")
# , 
#        subtitle="Starting trees: 3,163 bootstrapped RAxML gene trees
# Filtered: Trees with mean BS >= 40% or 50%
# Rearranged: Entire set rearranged at 40-70% BS cutoff
# All trees were rooted with Notung prior to reconcilliation")
XX

ggsave(plot = XX,filename = "Figure-3.pdf", path="figures", width = 7, height=7)
ggsave(plot = XX,filename = "Figure-3.png", path="figures", width = 7, height=7)

# plot duplications on tree

library("tidyverse")
library("ape")

rm(list=ls())

dup.files <- list.files(path = "output", pattern = "duplication.txt", recursive = TRUE, full.names = TRUE)
dup.files <- dup.files[grep(pattern = "obsolete", x = dup.files, perl = TRUE, invert = TRUE)]
nnames <- strsplit(dup.files, "/") %>% unlist(., recursive=FALSE) %>% matrix(., nrow=7*3+2, byrow=TRUE) %>% .[,2] %>% gsub("-", "_", .)
dup.tables <- lapply(dup.files, function(x) read.table(x, header=T, row.names = 1, na.strings = "NaN", sep = "\t"))
for (i in seq_along(nnames)) dup.tables[[i]]$Set <- nnames[i]
dup.tab.df <- bind_rows(dup.tables)

dup.tab.df3 <- dup.tab.df %>% 
  group_by(Set) %>% 
  #select(Set, starts_with("n", ignore.case = FALSE)) %>%
  mutate_all(function(x) x > 0) %>% 
  mutate(N=n()) %>% 
  ungroup %>% 
  group_by(Set, N) %>% 
  summarise_all(function(x) sum(x)) %>% 
  gather(Names, Count, -N, -Set) %>% 
  mutate(Percent=Count/N*100) %>% 
  dplyr::arrange(Names, .by_group=TRUE) %>% 
  mutate(Above15=ifelse(Percent >= 15, "yes", "no"))

ggplot(dup.tab.df3, aes(x=Names, y=Percent, group=Set, colour=Above15)) +
  geom_point() +
  facet_wrap("Set", ncol=1)

ST <- read.tree("output/notung-ST-all.newick")
plot(ST); nodelabels(ST$node.label)
sum(!ST$node.label %in% names(dup.tab.df))
sum(!ST$tip.label %in% names(dup.tab.df))
ST <- ladderize(ST)

pennST <- read.tree("output/notung-ST-penns.newick")
plot(pennST); nodelabels(pennST$node.label)
sum(!pennST$node.label %in% names(dup.tab.df))
sum(!pennST$tip.label %in% names(dup.tab.df))
pennST <- ladderize(pennST)

thalST <- read.tree("output/notung-ST-thals.newick")
plot(thalST); nodelabels(thalST$node.label)
sum(!thalST$node.label %in% names(dup.tab.df))
sum(!thalST$tip.label %in% names(dup.tab.df))
thalST <- ladderize(thalST)

# function to plot duplications on tree
plotTreeDups <- function(Tree, Dups, Tips=TRUE) {
  # Dups is a dataframe with columns Set, Names, Percent
  par(mar=c(0.5,0.5,1.5,0.5))
  plot.phylo(Tree, cex=.7, font=1, label.offset = 6)
  Title <- Dups["Set"] %>% unique
  title(main = Title, font.main=1, adj=0.05)
  Dups <- setNames(Dups$Percent, Dups$Names)
  Nsize <- Dups[Tree$node.label]
  Tsize <- Dups[Tree$tip.label]
  nodelabels(pch=21, cex=Nsize/20, bg=alpha("forestgreen", 0.8))
  nodelabels(text=round(Nsize, 0), frame="n", bg='white', adj=c(1.5,1.5), cex=.8, font=2)
  tiplabels(pch=21, cex=Tsize/20, bg=alpha("forestgreen", 0.8))
}

# test
# xx <- dup.tab.df3 %>% filter(Set=="cladeAll_bs50") %>% select(Names, Percent)
# plotTreeDups(Tree = ST, Dups = xx, Tips = TRUE)

# error because a data frame is not returned
# but plots are made and saved
pdf("figures/cladeAll_dups.pdf", onefile=TRUE)
dup.tab.df3 %>% filter(grepl("cladeAll", Set)) %>% group_by(Set) %>% 
  do({
    plotTreeDups(Tree=ST, Dups=.)
    })
dev.off()

pdf("figures/cladeThals_dups.pdf", onefile=TRUE)
dup.tab.df3 %>% filter(grepl("cladeThals", Set)) %>% group_by(Set) %>% 
  do({
    plotTreeDups(Tree=thalST, Dups=.)
  })
dev.off()

pdf("figures/cladePenns_dups.pdf", onefile=TRUE)
dup.tab.df3 %>% filter(grepl("cladePenns", Set)) %>% group_by(Set) %>% 
  do({
    plotTreeDups(Tree=pennST, Dups=.)
  })
dev.off()

### Ks bubbles
#library("ggjoy")
library("ape")
library("ggtree")
library("tidyverse")
library("viridis")

rm(list=ls())

## the tree
tr <- read.tree("original-data/absolute-time/absolute-time-rerooted-RT_38_72h.contree.tre")
tr <- ladderize(tr, right = TRUE)
# plot
plot.tr <- ggtree(tr=tr, ladderize = TRUE, right = FALSE) + 
  geom_tiplab() + 
  theme_tree2() +
  scale_x_continuous(breaks=rev(seq(250,0,-25)), labels=seq(250,0,-25))
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

Bubbles <- ggplot(data=Kstab, aes(group=Species2, y=Species2, x=Ks_value, colour=Pipeline, size=Paralog_pairs_under_peak/Total_paralog_pairs*100)) +
  geom_point(alpha=.8, stroke=.8, pch=24) + 
  #facet_wrap(~Species, ncol=1, scales="free_y", as.table = FALSE, strip.position = "left") +
  theme_minimal() +
  xlab("Ks") +
  ylab("") +
  scale_colour_viridis(discrete=TRUE) +
  guides(colour = guide_legend(override.aes = list(size=4))) +
  #scale_size_continuous(name="Paralog pairs\nunder peak", breaks = c(20, 600, 1000), labels = c(300, 600, 1000), range = c(0.5, 7)) +
  scale_size_continuous(name="Paralog pairs\nunder peak\n(% of total)", breaks = c(25, 45, 65), labels = c("25%", "45%", "65%"), range = c(1, 5)) +
  scale_x_continuous(breaks = seq(0,2,.5), labels = seq(0,2,.5), limits = c(0,2)) +
  theme(#legend.position="top",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.spacing.y = unit(.05, "cm"),
    axis.line.x=element_line(size=.35),
    axis.ticks.x=element_line(size=.35)
    # axis.text.y = element_blank()
  )
Bubbles

ggsave(Bubbles, filename = "paralogs-under-peak.pdf", path = "figures", width = 5, height=12)
YY <- cowplot::plot_grid(plot.tr, Bubbles, align = "hv", axis = "b", rel_widths = c(1,2))
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

#### Figure 5 summary

library("ape")
library("dplyr")

rm(list=ls())
Tree <- read.tree("original-data/absolute-time/absolute-time-rerooted-RT_38_72h.contree.tre")
Tree <- ladderize(Tree, right = FALSE)
DD <- data.frame(Treename=rev(Tree$tip.label),
                 Species=c("Triparma pacifica", "Corethron hystrix",
                           "Leptocylindrus danicus", "Proboscia alata", "Aulacoseira subarctica", "Stephanopyxis turris",
                           "Actinocyclus subtilis", "Rhizosolenia setigera", "Ditylum brightwellii", "Porosira pseudodelicatula",
                           "Conticribra weissflogii", "Cyclotella* nana", "Skeletonema marinoi", 
                           "Discostella pseudostelligera", "Thalassiosira oceanica", 
                           "Lampriscus shadboltianum", "Minutocellus polymorphus", "Odontella aurita",
                           "Eucampia antarctica", "Hemiaulus sinensis", "Chaetoceros affinis", "Chaetoceros sp.",
                           "Attheya septentrionalis", "Striatella unipunctata", "Asterionellopsis glacialis", "Talaroneis poseidonae",
                           "Staurosira sp.", "Diatoma tenuis", "Fragilaria sp.", "Thalassiothrix antarctica",
                           "Eunotia naegeli", "Cylindrotheca closterium", "Fragillariopsis* cylindrus", "Nitzschia sp.",
                           "Gyrosigma fasciola", "Phaeodactylum* tricornutum", "Craticula ambigua", "Sellaphora sp."), stringsAsFactors = FALSE)
DD$Genus <- strsplit(DD$Species, " ") %>% do.call(rbind, .) %>% .[,1]

#manually label edges
Edges <- rep(NA, nrow(Tree$edge))
ks <- setNames(c(73, 72, 66, 65, 57, 36, 35, 34, 33, 30, 18), 
               c("COHY", "LEDA", "RHSE", "ACSU", "THOC", "ATSE", "STUN", "TAPO", "ASGL", "DITE", "GYFA"))
ks <- ks[!(names(ks) %in% c("COHY", "LEDA", "RHSE", "ACSU", "ATSE", "STUN", "DITE"))]
ks <- c("GYFA"=69, "ASGL"=50, "TAPO"=51, "THOC"=43)
gc <- c("GYFA"=69, "ASGL.TAPO"=49, "ASGL"=50, "TAPO"=51, "ALL.BUT.STUN"=48, "MED.PEN"=16, "MED.PEN.MEL"=6, "THOC"=43)
rc <- c("MED.PEN.MEL"=6, "MED.PEN"=16, "PENN"=48, "THAL"=35)
gr <- c("THAL"=39, "PENNa"=48, "PENNb"=52, "MED.PEN"=16, "MED.PEN.MEL"=6)

Cols <- viridis::viridis(n=4,alpha = .8)
Nsp <- length(Tree$tip.label)  

# plot
plotFunc <- function(...) {
#  layout(matrix(1:2, ncol=2), widths=c(0.3,0.2))
  
  dpar <- par()
  par(mar=c(2,0.5,2,0.0), lend=2, lwd=.4)
  plot.phylo(Tree, edge.width = 3, edge.color = "grey60", font=1, label.offset=4, show.tip.label = TRUE)
  axis(side=1, line=-1, at=seq(250,0,-25), labels=seq(0,250,25), tck=-0.01, cex.axis=.8, mgp=c(3, 0.2, 0))
  mtext(side = 1, line=0.4, text = "Mya", cex=.9)
  # axis(side = 3, line=-1, at=seq(250,0,-25), labels=seq(0,250,25), tck=-0.01, cex.axis=.8, mgp=c(3, 0.2, 0))
  # mtext(side = 3, line=0.4, text = "Mya", cex=.9)
  edgelabels(edge=ks, pch=22, cex=2.5, frame="none", bg=Cols[1])
  edgelabels(edge=ks, text = "KS", cex=0.6, frame="none", bg=Cols[1], font=2, col="white")
  edgelabels(edge=gc, pch=22, cex=2.5, frame="none", bg=Cols[2], adj=15)
  edgelabels(edge=gc, text = "GC", cex=0.6, frame="none", bg=Cols[2], adj=-1, font=2, col="white")
  edgelabels(edge=rc, pch=22, cex=2.5, frame="none", bg=Cols[3], adj=0.5)
  edgelabels(edge=rc, text = "RS", cex=0.6, frame="none", bg=Cols[3], adj=0.5, font=2, col="black")
  edgelabels(edge=gr, pch=22, cex=2.5, frame="none", bg=Cols[4], adj=-14)
  edgelabels(edge=gr, text = "RM", cex=0.6, frame="none", bg=Cols[4], adj=1.95, font=2, col="black")
  legend(pch=22, pt.bg = Cols, x="topleft", bty="n", pt.cex=2.5, cex=.9,
         legend = c("Synon. distance", "Gene count", "Reconc. singly-labeled", "Reconc. multiply-labeled"))
#  text(x=rep(15,4), y = c(0.5,1.5,2.5,3.5)+2, labels = c("RM", "RS", "GC", "KS"), cex=.8, col=c("white", rep("black",3)), font = 2)
  
  # # labels
  # par(mar=c(2,0.0,2,0.0))
  # plot(y=1:Nsp, x=rep(0,Nsp), xaxt="n", yaxt="n", type = "n", bty="n", xlab="")
  # # axis(side=2, line=-0.5,at = 1:Nsp, labels=rep("", Nsp), las=2,
  # #      tcl=-0.2, mgp=c(3,0.5,0), cex.axis=.8, lwd=.2)
  # # axis(side=2, line=-0.5,at = 1:Nsp, labels=rep("", Nsp), las=2,
  # #      tcl=0.2, mgp=c(3,0.5,0), cex.axis=.8, lwd=.2)
  # text(y=1:Nsp, x=rep(-0.90,Nsp), labels = DD$Species, adj=0, cex=1.0)
  # par(lwd=0.6)
  # 
  # text(x=-0.9, y=0, "* Data from genomes", cex=.85, adj=0)
  # 
  # #  par(dpar)
}

pdf("figures/Figure-5.pdf", width=5, height=7)
plotFunc()
dev.off()

#####
# Figure 4

library("ape")
library("phangorn")

deepTreeA <- read.tree(text=c("(Triparma, (Corethron, (Leptocylindrus, (melosiroid, (coscinodiscoid, (multipolar, pennate))))));"))
deepTreeC <- read.tree(text=c("(Triparma, (Corethron, (Leptocylindrus, (melosiroid, (coscinodiscoid, (multipolar, pennate))))));"))
thalTree  <- read.tree(text=c("(Porosira, (Conticribra, (T_pseudonana, (Skeletonema, (Discostella, T_oceanica)))));"))
pennTree  <- read.tree(text=c("(multipolar_part,(Attheya, (Striatella, ((Asterionellopsis, Talaroneis), (araphid_rest, raphid)))));"))

plot(pennTree)
edgelabels()
deepEdgesA <- 6:12
deepEdgesC <- 10:12
thalEdges <- 6:10
pennEdges <- 4:10

deepPolyA <- 4:7
deepPolyC <- 6:7
thalPoly <- 4:6
pennPoly <- 4:7

deepH1A <- phangorn::mrca.phylo(x = deepTreeA, node = deepPolyA)
deepH2A <- c(2,4,3)
deepH2Af50 <- c(2,4,5)

deepH1C <- phangorn::mrca.phylo(x = deepTreeA, node = deepPolyC)
deepH2C <- c(2,4,5)
deepH2Cf50 <- c(2,4,5)

thalH1  <- phangorn::mrca.phylo(x = thalTree, node = thalPoly)
thalH2  <- c(2,4)
thalH2f50  <- c(2,4)

pennH1E <- phangorn::mrca.phylo(x = pennTree, node = pennPoly)
pennH2E <- c(2)
pennH2Ef50 <- c()

pennH1F <- phangorn::mrca.phylo(x = pennTree, node = pennPoly[3:4])
pennH2F <- c(2,3,5)
pennH2Ff50 <- c(5,7,6)

# plot function

Cols <- viridis::viridis(n=3, alpha = .7)

plotPhy <- function(Tree, h1, h2, h12, h22, Title="") {
  par(mar=rep(1,4))
  plot.phylo(Tree, font=1, label.offset=.1, type="p", cex=.8)
  xx <- rep(NA, Nnode(Tree))
  xx[h1-length(Tree$tip.label)] <- Cols[1]
  nodelabels(pch=21, bg=xx, col=xx, cex=2)
  edgelabels(pch=21, edge = h2, bg=Cols[2], cex=2)
  
  yy <- rep(NA, Nnode(Tree))
  yy[h12-length(Tree$tip.label)] <- alpha(Cols[1],.3)
  nodelabels(pch=21, bg=yy, col=yy, cex=2, adj=-.2)
  edgelabels(pch=21, edge = h22, bg=alpha(Cols[3],.3), cex=2, adj=-.2)
  
  title(Title, adj=0.0, font=1, cex.main=.85, line=-0.2)
}

pdf("figures/grampa-summary.pdf", width=8, height=2)
par(mfrow=c(1,5))
# unaltered and filter 50% bootstrap
plotPhy(deepTreeA, h1=deepH1A, h2=deepH2A, h12=deepH1A, h22=deepH2Af50, Title="Branch A")
legend(x=-1,y=7, pch=21, pt.bg=Cols, legend=c("h1", "h2 orig.", "h2 BS>50%"), bty="n", pt.cex=2)
plotPhy(deepTreeC, h1=deepH1C, h2=deepH2C, h12=deepH1C, h22=deepH2Cf50, Title="Branch C")
plotPhy(thalTree,  h1=thalH1,  h2=thalH2,  h12=thalH1,  h22=thalH2f50,  Title="Branch D'")
plotPhy(pennTree,  h1=pennH1E, h2=pennH2E, h12=pennH1E, h22=pennH2Ef50, Title="Branch E")
plotPhy(pennTree,  h1=pennH1F, h2=pennH2F, h12=pennH1F, h22=pennH2Ff50, Title="Branch F")
dev.off()

