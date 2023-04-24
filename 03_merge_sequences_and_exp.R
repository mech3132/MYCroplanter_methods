#!/bin/bash

### Merge sanger and experiments

dir.create("03_merge_sequences_and_exp")
library(ape)
library(tidyverse)

tre <- read.tree("02_process_sanger_illumina/tree_filtered_original.tre")
mtre <- read.tree("02_process_sanger_illumina/tree_filtered.tre")
dat <- read.delim("02_process_sanger_illumina/all_isolate_data.txt")
exp <- read.delim("2022-06-14-analyze_data_so_far/dat_clean.txt")

exp_singles_col1 <- exp %>% group_by(plate, isolate, experiment, Type) %>%
  summarise(root_stunt=mean(root_stunt, na.rm=TRUE)
            , whole_plant_sick=mean(whole_plant_sick, na.rm=TRUE)
            , growth_on_plant =mean(growth_on_plant, na.rm=TRUE)
            , Weight = mean(Weight, na.rm=TRUE)
            , growth_on_plate= mean(growth_on_plate, na.rm=TRUE)
            , k=mean(k, na.rm=TRUE)
            , r=mean(r, na.rm=TRUE)) %>%
  filter(Type%in%c("Single","MOCK")) %>% ungroup()
mock_col <- exp_singles_col1 %>% filter(Type=="MOCK") %>%
  select(plate, experiment, Weight) %>% rename(mockWeight=Weight)

exp_singles_col2 <- exp_singles_col1 %>%
  left_join(mock_col) %>%
  mutate(adjWeight= Weight-mockWeight) %>%
  group_by(isolate, Type) %>%
  summarise(root_stunt=mean(root_stunt, na.rm=TRUE)
            , whole_plant_sick=mean(whole_plant_sick, na.rm=TRUE)
            , growth_on_plant =mean(growth_on_plant, na.rm=TRUE)
            , adjWeight = mean(adjWeight, na.rm=TRUE)
            # , Weight_se = sd(Weight, na.rm=TRUE)/sqrt(n())
            , growth_on_plate= mean(growth_on_plate, na.rm=TRUE)
            , k=mean(k, na.rm=TRUE)
            , r=mean(r, na.rm=TRUE)
            , n=n()) %>%
  ungroup() 

# aveMock <- exp_singles_col2 %>% filter(isolate=="MOCK") %>% pull(Weight)
# exp_singles_col2_zeroed <- exp_singles_col2 %>% mutate(aveMock = aveMock) %>%
#   mutate(adjWeight=Weight-aveMock)

# Plot to see
exp_singles_col2 %>% 
  arrange(adjWeight) %>%
  mutate(isolate=factor(isolate, levels=unique(isolate))) %>%
  ggplot() + geom_bar(aes(x=isolate, y=adjWeight, fill=growth_on_plant), stat="identity") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

## Merge with seq dat
allDat <- exp_singles_col2 %>%
  rename(isolateID = isolate) %>% left_join(dat) %>%
  filter(Feature.ID %in% tre$tip.label)
trefilt <- keep.tip(tre, tip=allDat$Feature.ID)


# Get various colours
treetips <- trefilt$tip.label



NewLab <- pull(allDat[match(treetips, allDat$Feature.ID),"adjWeight"])


plot(trefilt)

###### For Sarzana

tree_filt_myc <- keep.tip(mtre, mtre$tip.label[grep("^AT",mtre$tip.label)])
# Get colors
isoNames <- pull(separate(data=data.frame(tipLabel=tree_filt_myc$tip.label),col=tipLabel, sep="_", into=c("isolateID")))
weightVals <- pull(allDat[match(isoNames, allDat$isolateID),"adjWeight"])
allDat %>% View()


colfunc <- colorRampPalette(c("red","grey", "blue"))
max10 <- round(max(abs(weightVals), na.rm = TRUE)*10)*2


colGrad <- data.frame(val=seq(-max10/2, max10/2), col= colfunc(max10+1))
colMapped <- colGrad[match(round((weightVals)*10), colGrad$val),"col"]
is.na(colMapped) <- "black"

png("03_merge_sequences_and_exp/tre_myc_withgrowth.png", height=10, width=8, unit="in", res=300)
plot(tree_filt_myc, tip.color = colMapped, cex=0.5)
dev.off()

######
