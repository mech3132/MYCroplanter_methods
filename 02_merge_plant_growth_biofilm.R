library(tidyverse)

dir.create("02_merged_plant_growth_biofilm")
######## Load data ######
plant <- read.delim("01b_merge_96HT/merged_plant_data.txt")
gc <- read.delim("01c_growthcurve_process_96HT/allgrowthfits_summary.txt")
cv <- read.delim("01d_biofilm_process/cv_data_summary.txt")

# All protect
allProtect <- plant %>% select(protect) %>% 
  filter(!protect %in% c("MOCK", NA)) %>%
  unique() %>% pull()

# All path
allPath <- plant %>% select(path) %>% 
  filter(!path %in% c("MOCK", NA)) %>%
  unique() %>% pull()

gc_cv <- full_join(gc, cv) %>%
  mutate(isolate = ifelse(isolate=="ATD3B1_1", "ATD3B1", isolate)) %>%
  filter()

gc_cv %>%
  ggplot() + 
  geom_point(aes(x=kmean, y=medianOD)) +
  geom_text(aes(x=kmean, y=medianOD, label=isolate))

# Get zeros/mock to subtract/adjust
plant_filt <- plant %>%
  filter(!is.na(path), !is.na(protect))
plant_gt <- plant %>%
  filter(experiment=="PolymyxaGT_2023-03-15", is.na(protect))
control_plants <- plant_filt %>%
  filter(Treatment=="MOCK-MOCK") %>%
  group_by(experiment, plate) %>%
  summarise(controlPixels = median(pixels)
            , controlGA = median(greenarea))

plant_adj <- full_join(plant_filt, control_plants) %>%
  mutate(adjPixels = pixels/controlPixels, adjGA = greenarea/controlGA)
write.table(plant_adj, file="02_merged_plant_growth_biofilm/plant_data.txt", row.names = FALSE, quote = FALSE, sep="\t")
write.table(plant_gt, file="02_merged_plant_growth_biofilm/plant_gt.txt", row.names = FALSE, quote = FALSE, sep="\t")

plant_summary <- plant_adj %>%
  filter(!is.na(adjPixels)) %>%
  group_by(Treatment, protect, path, ratio_PR_TO_PA, od_path, od_protect) %>%
  summarise(medianPix = median(adjPixels), sePix = sd(adjPixels)/n(), medianGA = median(adjGA), medianRawPixels = median(pixels), medianRawGA = median(greenarea)) %>%
  ungroup()
  
plant_summary_wgccv <- gc_cv %>%
  filter(isolate %in% allProtect) %>%
  rename_all(~paste0(.,"_protect")) %>%
  rename(protect = isolate_protect) %>%
  full_join(plant_summary, multiple = "all") %>%
  full_join(gc_cv %>%
              filter(isolate %in% allPath) %>%
              rename_all(~paste0(.,"_path")) %>%
              rename(path = isolate_path))

write.table(plant_summary_wgccv, file="02_merged_plant_growth_biofilm/merged_plant_gc_cv.txt", row.names = FALSE, quote = FALSE, sep="\t")















plant_summary_wgccv %>%
  ggplot() + 
  geom_point(aes(x=medianOD_protect, y=medianPix))


######### Effect of isolates ALONE ############
############ OD
plant_summary_wgccv %>%
  filter(path=="MOCK", protect!="ATD3B1") %>%
  arrange(medianPix) %>%
  mutate(protect = factor(protect, levels=unique(protect))) %>%
  ggplot() +
  geom_bar(aes(x=protect, y=medianPix, fill=medianOD_protect), stat='identity') +
  geom_errorbar(aes(x=protect, ymin=medianPix-sePix, ymax=medianPix+sePix))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

############ r
plant_summary_wgccv %>%
  filter(path=="MOCK", protect!="ATD3B1") %>%
  arrange(medianPix) %>%
  mutate(protect = factor(protect, levels=unique(protect))) %>%
  ggplot() +
  geom_bar(aes(x=protect, y=medianPix, fill=rmean_protect), stat='identity') +
  geom_errorbar(aes(x=protect, ymin=medianPix-sePix, ymax=medianPix+sePix))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

############ k
plant_summary_wgccv %>%
  filter(path=="MOCK", protect!="ATD3B1") %>%
  arrange(medianPix) %>%
  mutate(protect = factor(protect, levels=unique(protect))) %>%
  ggplot() +
  geom_bar(aes(x=protect, y=medianPix, fill=kmean_protect), stat='identity') +
  geom_errorbar(aes(x=protect, ymin=medianPix-sePix, ymax=medianPix+sePix))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))



######### Effect of isolates IN PAIRS ############
############ OD
plant_summary_wgccv %>%
  filter(path!="MOCK", protect!="ATD3B1", path!="N2C3") %>% 
  arrange(medianPix) %>%
  mutate(protect = factor(protect, levels=unique(protect))) %>%
  ggplot() +
  geom_bar(aes(x=protect, y=medianPix, fill=medianOD_protect), stat='identity') +
  geom_errorbar(aes(x=protect, ymin=medianPix-sePix, ymax=medianPix+sePix))+
  geom_hline(aes(yintercept=1), col="red")+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
    facet_grid(path~.)

############ r
plant_summary_wgccv %>%
  filter((path!="MOCK" ), protect!="ATD3B1", path!="N2C3") %>% 
  arrange(medianPix) %>%
  mutate(protect = factor(protect, levels=unique(protect))) %>%
  ggplot() +
  geom_bar(aes(x=protect, y=medianPix, fill=rmean_protect), stat='identity') +
  geom_errorbar(aes(x=protect, ymin=medianPix-sePix, ymax=medianPix+sePix))+
  geom_hline(aes(yintercept=1), col="red")+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  facet_grid(path~.)

############ k
plant_summary_wgccv %>%
  filter(path!="MOCK", protect!="ATD3B1", path!="N2C3") %>% 
  arrange(medianPix) %>%
  mutate(protect = factor(protect, levels=unique(protect))) %>%
  ggplot() +
  geom_bar(aes(x=protect, y=medianPix, fill=kmean_protect), stat='identity') +
  geom_errorbar(aes(x=protect, ymin=medianPix-sePix, ymax=medianPix+sePix))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  facet_grid(path~.)
