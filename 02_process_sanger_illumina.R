#!/bin/bash
library(ape)
library(tidyverse)
dir.create("02_process_sanger_illumina")

ill_tax <- read.delim("00_ILLUMINASEQUENCES/vsearch_silva138_16S_250-220/taxonomy.tsv") %>%
  rename_at(vars(everything()), ~paste0(., "_ill"))
ill_table <- read.delim("00_ILLUMINASEQUENCES/dada2_output_16S_250-220/feature-table.txt", skip=1)
san_dat <- read.csv("00_SANGERSEQUENCES/04_extract_projects/sequence_id_taxonomy_all.csv") %>%
  rename_at(vars(everything()), ~paste0(., "_san"))
# sansan_blast  <- read.delim("01_compare_illumina_sanger/blast_99_sanger_sanger/exported/blast6.tsv", header=FALSE) %>%
#   rename(Feature.ID_san = V1, Feature.ID_san2=V2)
# blast_results <-read.delim("01_compare_illumina_sanger/blast_97_ill_against_san/exported/blast6.tsv", header=FALSE) %>%
#   select(V1,V2,V3, V12) %>%
#   rename_at(vars(everything()), ~c("Feature.ID_ill", "Feature.ID_san","Identity","seqlength"))
blast_results <-read.delim("01_compare_illumina_sanger/blast_97/exported/blast6.tsv", header=FALSE) %>%
  select(V1,V2,V3, V12) %>%
  rename_at(vars(everything()), ~c("Feature.ID_san", "Feature.ID_ill","Identity","seqlength"))
tree <- read.tree("00_SANGERSEQUENCES/02_assigning_taxonomy/sepp_insertion/tree.nwk")


######## Merge illumina and sanger ###########

# Join all sequences
illsan_combined <- full_join(blast_results, ill_tax) %>%
  left_join(san_dat) 
# illsan_combined %>% View()

# Remove unassigned
illsan_combined_filt <- illsan_combined %>%
  select(Feature.ID_ill, Feature.ID_san, isolateID_san, Identity, seqlength, ScientificName_san, Taxon_san, Taxon_ill) %>%
  mutate(Taxon = ifelse(Taxon_san == "Unassigned" | is.na(Taxon_san), Taxon_ill, Taxon_san)) %>%
  filter(Identity>97) 

# List of unique Illumina sequences that have matches
rep_ill_seqs <- illsan_combined_filt %>% select(Feature.ID_ill) %>% distinct() %>%
  rename(X.OTU.ID = Feature.ID_ill) %>%
  mutate(Isolate=TRUE)

######### Taxa summary plot ##########
taxa_table_with_matches <- ill_table %>%
  select(X.OTU.ID, Mls1.16S_S264) %>%
  mutate(RelAbund = Mls1.16S_S264/sum(Mls1.16S_S264)) %>%
  left_join(rep_ill_seqs)
taxa_table_with_matches %>% filter(Isolate) %>% pull(RelAbund) %>% sum()

gg_rep_isolates <- taxa_table_with_matches %>%
  mutate(Isolate=ifelse(is.na(Isolate), FALSE, Isolate)) %>%
  ggplot() + 
  geom_bar(aes(x="Illumina\nSample", y=RelAbund, fill=Isolate), stat="identity") +
  xlab("") + ylab("Relative Abundance") +
  labs(fill=">97% isolate\nrepresentative") +
  scale_fill_manual(values=c("grey","green"))

ggsave(filename = "02_process_sanger_illumina/gg_rep_97_community.png", gg_rep_isolates,
       height=4, width=5)
taxa_table_with_matches %>% group_by(Isolate) %>% summarise(tot = sum(RelAbund))

#### Let's try collapsing by genus ####
ill_tax_adj <- ill_tax %>% 
  separate(Taxon_ill, sep ="; ",remove=FALSE, into=c("Kingdom_ill","Phylum_ill","Class_ill","Order_ill","Family_ill","Genus_ill","Species_ill"), fill="right") %>%
  rowwise() %>%
  mutate(Species = ifelse(length(grep("uncultured|sp[.]|metagenome|bacterium",Species_ill))>0, NA, Species_ill)) %>%
  mutate(Genus = ifelse(length(grep("Candidatus|uncultured|sp[.]|metagenome|bacterium|g__$",Genus_ill))>0, NA, Genus_ill)) %>%
  mutate(Family = ifelse(length(grep("Candidatus|uncultured|sp[.]|metagenome|bacterium|g__$",Family_ill))>0, NA, Family_ill)) %>%
  mutate(Order = ifelse(length(grep("Candidatus|uncultured|sp[.]|metagenome|bacterium|g__$",Order_ill))>0, NA, Order_ill)) %>%
  rename(X.OTU.ID = Feature.ID_ill) %>%
  ungroup() %>%
  select(X.OTU.ID, Species, Genus, Family,Order)

ill_table_withtax <- ill_table %>%
  select(X.OTU.ID, Mls1.16S_S264) %>%
  mutate(RelAbund = Mls1.16S_S264/sum(Mls1.16S_S264)) %>%
  left_join(ill_tax_adj)

## Samples that were identified to genus
gg_illumina <- ill_table_withtax %>%
  mutate(Genus_level=!is.na(Genus)) %>%
  ggplot() + 
  geom_bar(aes(x="Illumina\nSample", y=RelAbund, fill=Genus_level), stat="identity")+
  xlab("") + ylab("Relative Abundance") +
  labs(fill="Identified to\nGenera") +
  scale_fill_manual(values=c("grey","green"))
gg_illumina
ggsave(filename = "02_process_sanger_illumina/gg_illumina_to_genus.png", gg_illumina
       , height=4, width=5)

##### About 70% of Illuma sequences could be ID'ed to genus-level
isolate_genus <- san_dat %>% select(Genus_san) %>%
  distinct() %>% rename(Genus=Genus_san) %>%
  mutate(IsolateGenus=TRUE) 

ill_table_withtax %>%
  left_join(isolate_genus) %>%
  filter(IsolateGenus) %>% pull(RelAbund) %>% sum()

gg_genus_rep <- ill_table_withtax %>%
  left_join(isolate_genus) %>%
  # mutate(IsolateGenusRep = Genus %in% unique(san_dat$Genus_san)) %>%
  # mutate(Genus_level = !is.na(Genus)) %>%
  mutate(IsolateGenus = ifelse(is.na(IsolateGenus), FALSE, IsolateGenus)) %>%
  ggplot() + 
  geom_bar(aes(x="Illumina\nSample", y=RelAbund, fill=IsolateGenus), stat="identity") +
  xlab("") + ylab("Relative Abundance") +
  labs(fill="Genus rep\nin isolate \ncollection") + 
  scale_fill_manual(values=c("grey","green"))
gg_genus_rep
ggsave(filename = "02_process_sanger_illumina/gg_genus_rep.png", gg_genus_rep
       , height=4, width=5)

ill_table_withtax %>%
  group_by(Genus) %>%
  summarise(RelAbund = sum(RelAbund)) %>%
  filter(Genus=="g__Paenibacillus")

######### composition
ill_table_genus <- ill_table_withtax %>% group_by(Genus) %>%
  summarise(RelAbund = sum(RelAbund)) %>%
  mutate(Gr5 = !(RelAbund<0.01 | is.na(Genus)))
length(unique(ill_table_withtax$Genus))
high_genus <-ill_table_genus %>% filter(Gr5) %>%
  pull(Genus)
high_genus_colors <- c("maroon","blue","darkred","darkgreen","goldenrod","purple","brown","yellow","magenta","red","orange","pink","green")
names(high_genus_colors) <- high_genus
gg_comp <- ill_table_genus %>%
  arrange(-RelAbund) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus))) %>%
  ggplot() + geom_bar(aes(x="IlluminaSample", y=RelAbund, fill=Genus), stat='identity') +
  scale_fill_manual(values=high_genus_colors) +
  ylab("Relative Abundance") +xlab("")
gg_comp
ggsave(filename = "02_process_sanger_illumina/gg_genus_composition.png", gg_comp
       , height=4, width=5)


########## Looking at tree representation ############

san_dat


# Find duplicates
dat_good_qual <- san_dat %>% rowwise() %>%
  mutate(quality=ifelse(QualitySCore_san>=40& CRL_san>=500, 3, 
                        ifelse(QualitySCore_san>=25 | CRL_san>=500,2,1))) %>%
  ungroup() %>%
  filter(quality>=2)


dat_all_qual <- san_dat %>% rowwise() %>%
  mutate(quality=ifelse(QualitySCore_san>=40& CRL_san>=500, 3, 
                        ifelse(QualitySCore_san>=25 | CRL_san>=500,2,1))) %>%
  ungroup()

dupIsos <- dat_all_qual  %>%
  mutate(dup = duplicated(isolateID_san)) %>%
  filter(dup) %>% pull(isolateID_san)


duplegend <- dat_all_qual %>% select(Project_san, Feature.ID_san, isolateID_san, ScientificName_san, quality) %>%
  # filter(isolateID %in% dupIsos) %>%
  group_by(isolateID_san) %>% mutate(dupN=rank(isolateID_san, ties.method = "random")) %>% ungroup() 
dat_all_singles <- duplegend %>% select(-Feature.ID_san) %>%
  pivot_wider(names_from=dupN, values_from=c(ScientificName_san,quality)) %>%
  rowwise() %>%
  mutate(whichM = which.max(c(quality_1, quality_2,quality_3))) %>% 
  mutate(ScientificName = get(paste0("ScientificName_san_",whichM)), quality=get(paste0("quality_",whichM)), dupN=whichM) %>% 
  select(isolateID_san, ScientificName,quality, dupN)
# Add back feature ids
dat_all_singles_comp <- duplegend %>% select(Project_san, Feature.ID_san, isolateID_san, dupN) %>%
  right_join(dat_all_singles) %>%
  filter(Feature.ID_san %in% tree$tip.label)

######## Tree
# Remove faile dtree insertions
tree_filt_original <- keep.tip(tree, dat_all_singles_comp$Feature.ID_san)
tree_filt <- keep.tip(tree, dat_all_singles_comp$Feature.ID_san)

projCol <- ifelse(pull(dat_all_singles_comp[match(tree_filt$tip.label, dat_all_singles_comp$Feature.ID_san),"Project_san"])=="MYC"
                  ,"red","blue")

isoNames <- pull(dat_all_singles_comp[match(tree_filt$tip.label, dat_all_singles_comp$Feature.ID_san),"isolateID_san"])
isoIDs <- pull(dat_all_singles_comp[match(tree_filt$tip.label, dat_all_singles_comp$Feature.ID_san),"ScientificName"])
newNames <- paste0(isoNames," ",isoIDs)
tree_filt$tip.label <- newNames


png("02_process_sanger_illumina/all_iso_tree.png", height=12, width=6, units="in", res=300)
plot(tree_filt, cex=0.5, tip.color = projCol)
dev.off()

write.tree(tree_filt, file = "02_process_sanger_illumina/tree_filtered.tre")
write.tree(tree_filt_original, file = "02_process_sanger_illumina/tree_filtered_original.tre")

allIsos <- right_join(dat_all_qual, dat_all_singles_comp) 
colnames(allIsos) <- gsub("_san","",colnames(allIsos))

write.table(allIsos, file="02_process_sanger_illumina/all_isolate_data.txt", row.names = FALSE
            , quote = FALSE, sep="\t")
