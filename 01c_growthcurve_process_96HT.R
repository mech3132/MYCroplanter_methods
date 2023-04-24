library(tidyverse)
library(growthcurver)

dir.create("01c_growthcurve_process_96HT")
#### Clean individual growth curve data ########
# Too large to do all at once, so I'm calcuating growth curves separately for each #
# BATCH LIST 
batchlist <- c(1,2,3,4)
##### 2022-04-28 #######
gc1 <- read.csv("2022-04-28_rhizo_iso_screen/00_data/platereader/2022-04-28_rhizo_iso_growthcurves.csv", na.strings = c("","NA","na")) %>%
  separate(time, into = c("d","other"), sep="[.]", fill = "left") %>%
  separate(other, into = c("h","m","s"), sep=":") %>%
  mutate(d = ifelse(is.na(d),0,as.numeric(d))) %>%
  mutate(time = d*24 + as.numeric(h) + as.numeric(m)/60 ) %>%
  select(experiment, time, temperature, A1:H12)
meta1 <- read.csv("2022-04-28_rhizo_iso_screen/00_meta/2022-04-28_iso_layout_meta.csv", na.strings = c("","NA","na"))
##### 2022-05-11 #######
gc2 <- read.csv("2022-05-11_growthcurve_aileen/00_dat/2022-05-11_plate_od.csv", na.strings = c("","NA","na"))%>%
  separate(time, into = c("d","other"), sep="[.]", fill = "left") %>%
  separate(other, into = c("h","m","s"), sep=":") %>%
  mutate(d = ifelse(is.na(d),0,as.numeric(d))) %>%
  mutate(time = d*24 + as.numeric(h) + as.numeric(m)/60 ) %>%
  select(experiment, time, temperature, A1:H12)
meta2 <- read.csv("2022-05-11_growthcurve_aileen/00_metadata/2022-05-11_plate_isos.csv", na.strings = c("","NA","na"))

##### 2022-05-24 #######
gc3 <- read.csv("2022-05-24_growthcurve/00_dat/growth_curve_isos_2022-05-25_format.csv")%>%
  separate(time, into = c("d","other"), sep="[.]", fill = "left") %>%
  separate(other, into = c("h","m","s"), sep=":") %>%
  mutate(d = ifelse(is.na(d),0,as.numeric(d))) %>%
  mutate(time = d*24 + as.numeric(h) + as.numeric(m)/60 ) %>%
  select(experiment, time, temperature, A1:H12)
meta3 <- read.csv("2022-05-24_growthcurve/00_meta/2022-05-25_plate_layout.csv")

##### 2023-03-24 #######
gc4 <- read.csv("2023-03-24_growthcurves/dat/2023-03-24_growthcurves.csv")%>%
  separate(time, into = c("d","other"), sep="[.]", fill = "left") %>%
  separate(other, into = c("h","m","s"), sep=":") %>%
  mutate(d = ifelse(is.na(d),0,as.numeric(d))) %>%
  mutate(time = d*24 + as.numeric(h) + as.numeric(m)/60 ) %>%
  select(experiment, time, temperature, A1:H12)
meta4 <- read.csv("2023-03-24_growthcurves/meta/meta_growthcurves_2023-03-24.csv")

####  Process all to standardize input info ####
colnames_for_gc <- c("experiment","time","temperature","well","rawod")
colnames_for_dat <- c("plate","isolate","well","experiment")

gc1_long <- gc1 %>%
  pivot_longer(-c(experiment, time, temperature), names_to = "well", values_to="rawod") %>%
  select(all_of(colnames_for_gc))
gc2_long <- gc2 %>%
  pivot_longer(-c(experiment, time, temperature), names_to = "well", values_to="rawod") %>%
  select(all_of(colnames_for_gc))
gc3_long <- gc3 %>%
  pivot_longer(-c(experiment, time, temperature), names_to = "well", values_to="rawod") %>%
  select(all_of(colnames_for_gc))
gc4_long <- gc4 %>%
  pivot_longer(-c(experiment, time, temperature), names_to = "well", values_to="rawod") %>%
  select(all_of(colnames_for_gc))

alldat1 <- meta1 %>% as_tibble() %>% rename(well=well_od) %>%
  mutate(plate=1) %>%
  select(all_of(colnames_for_dat)) %>%
  full_join(gc1_long, multiple = "all") %>%
  filter(!is.na(well))
alldat2 <- meta2 %>% as_tibble() %>% 
  select(all_of(colnames_for_dat)) %>%
  full_join(gc2_long, multiple="all")
alldat3 <- meta3 %>% as_tibble() %>%
  unite("row","col", col=well, sep="") %>%
  mutate(plate=1) %>%
  select(all_of(colnames_for_dat)) %>%
  full_join(gc3_long, multiple="all")
alldat4 <- meta4 %>% as_tibble() %>%
  unite("row","col", col=well, sep="") %>%
  select(all_of(colnames_for_dat)) %>%
  full_join(gc4_long, multiple="all") 

####### Process each dataset individually ###########
allFits <- data.frame()
for ( b in batchlist ) {
  tempdat <- get(paste0("alldat",b))
  tempmat <- tempdat %>%
    unite(isolate, well, plate,experiment, col="uniquewell", sep=":") %>%
    select(uniquewell, rawod, time) %>%
    pivot_wider(names_from=uniquewell, values_from=rawod)
  
  tempfit <- SummarizeGrowthByPlate(tempmat)
  allFits <- rbind(tempfit, allFits)
  
}
#edit result table to get isolateIDs
allFits_edit <- allFits %>%
  separate(sample, into=c("isolate","well","plate","experiment"), sep=":") 



write.table(allFits_edit, file="01c_growthcurve_process_96HT/allgrowthfits.txt", row.names = FALSE, quote = FALSE, sep="\t")

# Plot growth rates
allFits_edit %>%
  filter(note=="", r<6, k<10) %>% # remove questionable data
  ggplot(aes(x=isolate, y=r)) +
  geom_boxplot()+
  geom_jitter(height=0, width=0.25, col="yellow") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

# Plot growth rates by carrying capacity
allFits_summary <- allFits_edit %>%
  filter(note=="", r<6, k<10) %>% # remove questionable data
  group_by(isolate) %>%
  summarise(rmean=mean(r), kmean=mean(k), rse = sd(r)/n(), kse = sd(k)/n())
allFits_summary %>%
  ggplot(aes(x=rmean, y=kmean)) +
  geom_point() + 
  geom_errorbar(aes(ymin=kmean-kse, ymax=kmean+kse))+
  geom_errorbarh(aes(xmin=rmean-rse, xmax=rmean+rse))

write.table(allFits_summary, file="01c_growthcurve_process_96HT/allgrowthfits_summary.txt", row.names = FALSE, quote = FALSE, sep="\t")

control_data <- full_join(alldat1, alldat2) %>%
  full_join(alldat3) %>%
  full_join(alldat4) %>%
  filter(isolate =="CONT"| isolate=="BLANK",!is.na(well))

control_data %>%
  ggplot() + 
  geom_point(aes(x=time, y=rawod, col=experiment)) +
  facet_wrap(.~isolate)

# 
# ## Cobine; edit time ###### TOO BIG, should process each growth curve separately
alldat <- full_join(alldat1, alldat2) %>%
  full_join(alldat3) %>%
  full_join(alldat4) %>%
  filter(isolate_raw !="CONT", !is.na(well)) %>%
#   unite(well, experiment, col=uniquewell,  sep = "_",remove=FALSE) %>%
#   # filter(!is.na(time), !is.na(isolate)) %>%
#   separate(time, into=c("dd","other"),sep="[.]", remove=FALSE, fill="left") %>%
#   mutate(dd=ifelse(is.na(dd), 0,dd)) %>%
#   separate(other, into=c("hh","mm","ss"), sep=":") %>%
#   mutate(time_h = as.numeric(dd)*24 + as.numeric(hh) + as.numeric(mm)/60 + as.numeric(ss)/360)
# # Check on growth of controls; remove H12 because it is weird
# alldat %>%
#   filter(isolate%in%c("BLANK")) %>%
#   ggplot() + geom_line(aes(x=time_h, y=rawod, group=uniquewell, col=well))
#  
# # Get blanks to adjust out
# blankods <- alldat %>% 
#   filter(isolate == "BLANK", isolate=="CONT",well !="H12") %>% 
#   group_by(experiment, time_h) %>%
#   summarize(blankOD = median(rawod)) %>% ungroup() 
# 
# alldat_adj <- alldat %>% full_join(blankods) %>%
#   filter(!is.na(rawod)) %>%
#   select(isolate, uniquewell, experiment, time_h, rawod, blankOD) %>%
#   mutate(ODadj = rawod-blankOD) 
# 
# alldat_adj %>%
#   ggplot() +geom_line(aes(x=time_h, y=ODadj, group=uniquewell, col=experiment)) +
#   facet_wrap(.~isolate, nrow=10, scales = "free_y")
