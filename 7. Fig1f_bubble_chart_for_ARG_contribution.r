
list.files()

library(data.table)
library(dplyr)

varMapping_df <- fread("var_mapping.txt")
meta <- fread("meta.txt")

# cutoff -------------------------

co_taxa.relAbund <- 1e-4

co_numARG.ofTaxa <- 2
co_avgContribution <- 0.01

Taxa.remove <- c("Escherichia coli")

# first identify minor otus -------------

load("1.otu_info.RData")
taxon.data <- 
  fread("kraken_taxa.txt", data.table = F) %>%
  mutate(species = sapply(`#NAME`, function(x)kraken_otu.info$Species[rownames(kraken_otu.info)==x])) %>%
  select(-`#NAME`)

taxon.relAbund <- 
  taxon.data %>% 
  reshape2::melt(id.var = "species", variable.name="sample", value.name = "count") %>%
  mutate(sp.Grp = sapply(sample, function(x) meta$Group[meta$NAME == x])) %>%
  group_by(sp.Grp, species) %>%
  summarise(counts = sum(count)) %>% mutate(relAbund = counts/sum(counts)) %>% 
  arrange(desc(relAbund)) %>% arrange(sp.Grp)

Taxa.keep <- taxon.relAbund$species[taxon.relAbund$relAbund > co_taxa.relAbund & taxon.relAbund$sp.Grp == "C"]


# differentiated ARG   -------------
diffARGs.dat <- 
  fread("NEU_EOS_cor_mod.txt") %>%
 # filter(P < 0.1) %>%
 # mutate(Group=sapply(log2FC, function(x)if(x<0)"EOS_enriched"else"NEU_enriched"))%>%
  mutate(ARG = sapply(NAME, function(x)varMapping_df$Gene[varMapping_df$Var==x]))



load("2.integrated_GuangZ.ShenZ_ARGsOAP.Taxon.RData")

# calculate alpha diversity of each ARG from species contribution--------------
library(vegan)
ARG.taxa.matrix <- ARG_species_df %>% 
  filter(!is.na(taxa)) %>%  # remove unidentified 
  mutate(sp.Grp = sapply(sample_name, function(x)meta$Group[meta$NAME == x])) %>%
  filter(sp.Grp == "C") %>%
  # filter(ARG %in% diffARGs.dat$ARG) %>%
  group_by(ARG, taxa) %>%
  summarise(readCount = sum(readCount)) %>%
  reshape2::dcast(ARG~taxa, value.var = "readCount") %>%
  tibble::column_to_rownames("ARG")

ARG.taxa.matrix[is.na(ARG.taxa.matrix)] <- 0

diversity_df1 <- cbind.data.frame(shannon = diversity(ARG.taxa.matrix, index = 'shannon'),
                                  invsimp = diversity(ARG.taxa.matrix, index = 'invsimpson'),
                                  simp = diversity(ARG.taxa.matrix, index = 'simpson'),
                                  specnumber = specnumber(ARG.taxa.matrix))

write.table(diversity_df1, file = "diversity.dat.txt", sep = "\t", quote = F, row.names = T)

# calculate alpha diversity of each ARG from major species contribution--------------
ARG.taxa.matrix <- ARG_species_df %>% 
  filter(!is.na(taxa)) %>%  # remove unidentified 
  mutate(sp.Grp = sapply(sample_name, function(x)meta$Group[meta$NAME == x])) %>%
  filter(sp.Grp == "C") %>%
  filter(taxa %in% Taxa.keep) %>%
  filter(!taxa %in% Taxa.remove) %>%
  #filter(ARG %in% diffARGs.dat$ARG) %>%
  group_by(ARG, taxa) %>%
  summarise(readCount = sum(readCount)) %>%
  reshape2::dcast(ARG~taxa, value.var = "readCount") %>%
  tibble::column_to_rownames("ARG")

ARG.taxa.matrix[is.na(ARG.taxa.matrix)] <- 0

diversity_df2 <- cbind.data.frame(shannon = diversity(ARG.taxa.matrix, index = 'shannon'),
                                  invsimp = diversity(ARG.taxa.matrix, index = 'invsimpson'),
                                  simp = diversity(ARG.taxa.matrix, index = 'simpson'),
                                  specnumber = specnumber(ARG.taxa.matrix))

write.table(diversity_df2, file = "diversity_fromMajorSpecies.dat.txt", sep = "\t", quote = F, row.names = T)

# contribution of species to each ARG  -------------
taxComp.inARG_dat <- ARG_species_df %>% 
  filter(!is.na(taxa)) %>%
  filter(taxa %in% Taxa.keep) %>%
  filter(!taxa %in% Taxa.remove) %>%
  filter(ARG %in% diffARGs.dat$ARG) %>%
  mutate(sp.Grp = sapply(sample_name, function(x)meta$Group[meta$NAME == x])) %>%
  filter(sp.Grp == "C") %>%
  group_by(ARG, taxa) %>% 
  summarise(value = sum(rpkm)) %>%
  mutate(Contribution = value/sum(value))

# keep only species with contribution to at least co_numARG.ofTaxa ARGs, and a average contribution > co_avgContribution
taxa.stats <- 
  taxComp.inARG_dat %>% group_by(taxa) %>%
  summarise(numARG = n(), avgContribution = mean(Contribution))


plotDat <- taxComp.inARG_dat %>%
  filter(taxa %in% taxa.stats$taxa[taxa.stats$numARG > co_numARG.ofTaxa &
                                     taxa.stats$avgContribution > co_avgContribution]) %>%
  mutate(phylum = sapply(taxa,
                         function(x) kraken_otu.info$Phylum[kraken_otu.info$Species == x])) %>%
  mutate(ARGtype = sapply(ARG, function(x) varMapping_df$Group[varMapping_df$Gene == x]))


# ba ba tu -------------

library(ggplot2)

# color by phyla


# arrange species by phyla
taxa_order <- plotDat %>% as.data.frame() %>% select(taxa, phylum) %>% unique() %>% arrange(taxa) %>% arrange(phylum)
plotDat$taxa <- factor(plotDat$taxa, levels = taxa_order$taxa)

# arrange ARG by ARG type
#ARG_order <- plotDat %>% as.data.frame() %>% select(ARG, ARGtype) %>% unique() %>% arrange(ARG) %>% arrange(ARGtype)

plotDat$ARG <- factor(plotDat$ARG, levels = diffARGs.dat$ARG)

# arrange ARG by log2FC
#plotDat$log2FC = sapply(plotDat$ARG,
#                        function(x) diffARGs.dat$log2FC[diffARGs.dat$ARG == x])
#ARG_order <- plotDat %>% as.data.frame() %>% select(ARG, log2FC) %>% unique() %>% arrange(log2FC)
#plotDat$ARG <- factor(plotDat$ARG, levels = ARG_order$ARG)

p<- ggplot2::ggplot(data = plotDat) +
  geom_point(aes(x=ARG,y=taxa, size=Contribution, color=phylum), shape =21) +
  theme_bw() +  theme(axis.text.x  = element_text(angle = 90), 
                      panel.grid = element_blank())
p
ggsave(p, filename = "testBaBaPlot.pdf", device = "pdf", width = 10, height = 10 )
