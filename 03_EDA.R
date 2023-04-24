library(tidyverse)

dir.create("03_EDA")

#### load data ######
plant <- read.delim("02_merged_plant_growth_biofilm/plant_data.txt") %>%
  mutate(ScientificName = ifelse(protect=="MOCK", "",ScientificName)) %>%
  unite(protect, ScientificName, col="IsolateName", remove=FALSE) 
plant_gt <- read.delim("02_merged_plant_growth_biofilm/plant_gt.txt")
merged <- read.delim("02_merged_plant_growth_biofilm/merged_plant_gc_cv.txt")

## Visualize mocks ##
plant %>%
  filter(Treatment=="MOCK-MOCK") %>%
  ggplot() +
  geom_point(aes(x=experiment, y=greenarea, col=factor(plate)))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ylab("Green leaf area (sqrt pixels)") +
  xlab("Experiment") +labs(col="Plate")

### Binomially distributed pixels #####

## Bimodal
gg_bin <- plant %>%
  ggplot() +
  geom_histogram(aes(x=greenarea))+
  facet_wrap(.~experiment) +
  xlab("Experiment")

gg_bin_clustered <- plant %>%
  ggplot() +
  geom_histogram(aes(x=greenarea, group=clusterDA, fill=clusterDA), show.legend = FALSE)+
  facet_wrap(.~experiment)+
  xlab("Experiment")
ggsave("03_EDA/gg_histogram_GA.png", gg_bin, height=4, width=5)
ggsave("03_EDA/gg_histogram_GA_clustered.png", gg_bin_clustered, height=4, width=5)

#####  Verifying N2C3/WCS365 output  ###########
exp_to_keep_n2c3 <- plant %>%
  filter(path=="N2C3" | protect == "WCS365") %>%
  pull(experiment) %>% unique()
plant_n2c3wcs <- plant %>%
  filter(experiment %in% exp_to_keep_n2c3, path %in% c("MOCK","N2C3"), protect %in% c("MOCK","WCS365"))
gg_basicnsc3wcs <- plant_n2c3wcs %>%
  filter(plant_age_inoc==5) %>%
  ggplot() +
  geom_boxplot(aes(x=Treatment, y=greenarea)) +
  geom_jitter(aes(x=Treatment, y=greenarea), col="green")+
  # facet_wrap(.~plant_age_inoc)+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ylab("Green leaf area (sqrt pixels)")
gg_basicnsc3wcs
ggsave("03_EDA/gg_basicnsc3wcs.png", gg_basicnsc3wcs, height=4, width=6)


#####  Ratio of N2C3 and WCS365 determines dead/alive  ###########

gg_ratios_ga <- plant %>%
  filter(experiment=="Grad_2023-02-01") %>%
  # filter(od_path>0) %>%
  ggplot(aes(x=(od_path), y=greenarea)) +
  geom_point(height=0.1, width=0)+
  # geom_smooth(method="glm",method.args = list(family = "binomial"))+
  geom_smooth(method="lm")+
  geom_hline(aes(yintercept=0.5), col="red", lty=2)+
  facet_grid(.~od_protect) +
  xlab("OD of pathogen (N2C3)") +
  ylab("Green leaf area (sqrt pixels)") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) 
gg_ratios_ga
ggsave("03_EDA/gg_N2C3_WCS_ratios_GA.png", gg_ratios_ga, height=5, width=10)


gg_ratios <- plant %>%
  # filter(od_path>0) %>%
  filter(experiment=="Grad_2023-02-01") %>%
  ggplot(aes(x=(od_path), y=clusterDA)) +
  geom_jitter(height=0.1, width=0)+
  geom_smooth(method="glm",method.args = list(family = "binomial"))+
  geom_hline(aes(yintercept=0.5), col="red", lty=2)+
  facet_grid(.~od_protect) +
  xlab("OD of pathogen (N2C3)") +
  ylab("Dead or alive (points jittered)") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) 
gg_ratios
ggsave("03_EDA/gg_N2C3_WCS_ratios.png", gg_ratios, height=5, width=10)

### Creating heatmap
gg_survivalheatmap <- plant %>%
  # filter(od_path>0) %>%
  filter(experiment=="Grad_2023-02-01") %>%
  group_by(od_path, od_protect) %>%
  summarise(p_infect = sum(clusterDA)/n()) %>%
  ggplot() + 
  geom_contour_filled(aes(x=od_path, y=od_protect, z=p_infect)) +
  labs(fill="Probability survival") +
  xlab("Pathogen OD")+ylab("Protective OD") +
  xlim(0.005, 0.01)
gg_survivalheatmap
ggsave(filename = "03_EDA/gg_survivalheatmap.png", gg_survivalheatmap, width=5, height=3)

###### summarising outcomes #####
# Get ordered list of protectives
allProtect_ordered <- plant %>% 
  filter(path=="MOCK") %>%
  group_by(protect) %>%
  summarise(med = median(greenarea)) %>%
  arrange(-med) %>%
  pull(protect) %>% unique()
allProtectname_ordered  <- plant %>% 
  filter(path=="MOCK") %>%
  group_by(IsolateName) %>%
  summarise(med = median(greenarea)) %>%
  arrange(-med) %>%
  pull(IsolateName) %>% unique()

plant_adj <- plant %>%
  filter(path %in% c("MOCK","ATR2A1","ATY2A2")) %>%
  mutate(path = factor(path, levels=c("MOCK","ATR2A1","ATY2A2"))
         , protect=factor(protect, levels=allProtect_ordered)
         , IsolateName=factor(IsolateName, levels=allProtectname_ordered))
### Pathogen effects ###
gg_patheffects_5 <- plant_adj %>%
  filter(plant_age_inoc==5) %>%
  filter(protect=="MOCK") %>%
  ggplot(aes(x=path, y=greenarea)) + 
  geom_boxplot()+
  geom_jitter(aes(col=path), height=0, width=0.25, show.legend = FALSE) +
  ylab("Median green leaf area\n(sqrt pixels)") + xlab("Pathogen")+
  scale_color_manual(values=c("blue","red","orange"))+
  labs("Treatment")
gg_patheffects_5
ggsave("03_EDA/gg_path_effects_day5.png", gg_patheffects_5, width=5, height=4)

gg_patheffects_allage <- plant_adj %>%
  filter(plant_age_inoc!=6) %>%
  filter(protect=="MOCK") %>%
  ggplot(aes(x=path, y=greenarea)) + 
  geom_boxplot()+
  geom_jitter(aes(col=path), height=0, width=0.25, show.legend = FALSE) +
  ylab("Median green leaf area\n(sqrt pixels)") + xlab("Pathogen")+
  scale_color_manual(values=c("blue","red","orange"))+
  labs("Treatment")+
  facet_grid(.~plant_age_inoc)+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
gg_patheffects_allage
ggsave("03_EDA/gg_path_effects_all_ages.png", gg_patheffects_allage, width=5, height=3)

gg_patheffects <- plant_adj %>%
  # filter(plant_age_inoc==5) %>%
  filter(protect=="MOCK") %>%
  ggplot(aes(x=path, y=greenarea)) + 
  geom_boxplot()+
  geom_jitter(aes(col=path), height=0, width=0.25, show.legend = FALSE) +
  ylab("Median green leaf area") + xlab("Pathogen")+
  scale_color_manual(values=c("blue","red","orange"))+
  labs("Treatment")
gg_patheffects
ggsave("03_EDA/gg_path_effects_allages_pooled.png", gg_patheffects, width=5, height=4)

##### Individual effects ########
plant_adj %>%
  filter(path=="MOCK") %>%
  group_by(protect) %>%
  summarise(Alive = sum(clusterDA==1)/n(),
            Dead = sum(clusterDA==0)/n()) %>%
  arrange(-Alive) %>%
  mutate(protect = factor(protect, levels=unique(protect))) %>%
  pivot_longer(-protect, names_to="DA", values_to="prop") %>%
  ggplot() +
  geom_bar(aes(x=protect, y=prop, group=DA, fill=DA), stat="identity")+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) 

plant_adj %>%
  filter(path=="MOCK") %>%
  ggplot() +
  geom_bar(aes(x=protect, group=clusterDA, fill=clusterDA))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) 

######## Difference in effect
plant_adj %>%
  # group_by(protect, path) %>%
  # summarise(median_GA = median(greenarea), se_GA = sd(greenarea)/n()) %>%
  # ungroup() %>%
  ggplot(aes(x=protect, y=greenarea, col=path)) +
  geom_boxplot()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  facet_grid(plant_age_inoc ~.)
  

gg_pairedeffects <- plant_adj %>%
  # filter(plant_age_inoc==5) %>%
  unite(protect, ScientificName, col="RowLabs", remove=FALSE) %>%
  group_by(protect, RowLabs, path) %>%
  summarise(median_GA = median(greenarea, na.rm=TRUE), se_GA = sd(greenarea, na.rm=TRUE)/n()) %>%
  ungroup() %>%
  ggplot(aes(x=protect, y=median_GA, col=path)) +
  geom_point()+
  geom_errorbar(aes(ymin=median_GA-se_GA, ymax=median_GA+se_GA)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ylab("Median green leaf area") + xlab("Candidate protector")+
  labs(col="Pathogen")+
  scale_color_manual(values=c("blue","red","orange"))
gg_pairedeffects
ggsave("03_EDA/gg_paired_effects.png", gg_pairedeffects, width=10, height=4)
  

gg_pairedeffects2 <- plant_adj %>%
  group_by(protect, IsolateName, path) %>%
  summarise(median_GA = median(greenarea, na.rm=TRUE), se_GA = sd(greenarea, na.rm=TRUE)/n()) %>%
  ungroup() %>%
  ggplot(aes(x=IsolateName, y=median_GA, col=path)) +
  geom_point()+
  geom_errorbar(aes(ymin=median_GA-se_GA, ymax=median_GA+se_GA)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ylab("Median green leaf area") + xlab("Candidate protector")+
  labs(col="Pathogen")+
  scale_color_manual(values=c("blue","red","orange"))

gg_pairedeffects2
ggsave("03_EDA/gg_paired_effects_names.png", gg_pairedeffects2, width=10, height=4)
length(allProtect_ordered)

## Just ones with multiple ages

treat_to_keep <- plant_adj %>%
  group_by(path, protect,Treatment) %>%
  summarise(maxAge = max(plant_age_inoc), minAge = min(plant_age_inoc)) %>%
  filter(maxAge==7 & minAge==5) %>%
  ungroup() %>%
  # filter(protect!="MOCK" | (protect=="MOCK" & path=="MOCK")) %>%
  pull(Treatment)

gg_ageffect <- plant_adj %>%
  filter(Treatment %in% treat_to_keep,plant_age_inoc!=6 ) %>%
  # filter(path=="MOCK"| protect=="MOCK") %>%
  mutate(iso = ifelse(path=="MOCK",paste0(protect), paste0(path))) %>%
  group_by(iso) %>%
  mutate(min_greenarea = min(greenarea)) %>%
  ungroup() %>%
  arrange(-greenarea) %>%
  mutate(iso = factor(iso, levels=unique(c("MOCK",iso)))) %>%
  ggplot(aes(x=iso, y=greenarea, fill=as.factor(plant_age_inoc))) +
  geom_boxplot()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_fill_manual(values=c("green","darkgreen"))+
  labs(fill="Plant age\nat inoculation\n(days)") +
  ylab("Green leaf area (sqrt pixels)") +
  xlab("Isolate (monoculture)")
gg_ageffect
ggsave("03_EDA/gg_age_effect.png", gg_ageffect, width=6, height=4)

# allProtect
plant_adj %>%
  group_by(protect, path) %>%
  mutate(propAlive = sum(clusterDA)/n()) %>%
  arrange(propAlive) %>%
  mutate(protect = factor(protect, levels=allProtect_ordered))%>%
  ggplot(aes(x=protect, y=propAlive, fill=path, group=path)) +
  geom_bar(stat="identity", position = "dodge")+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) 

###### Investigate trends ########
#Isolatesofinterest
isointerest <- c("ATY4B2","ATL3B2","ATD4A1","ATD3B6","ATR3A1")
# make filtered set
merged_isointerest <- merged %>%
  filter(protect%in%isointerest) %>%
  select(protect, rmean_protect, kmean_protect, rse_protect, kse_protect, medianOD_protect) %>%
  distinct()

gg_ecology_protect<- merged %>%
  filter(path%in%c("ATR2A1","ATY2A2")) %>%
  ggplot(aes(x=rmean_protect, y=kmean_protect)) + 
  geom_point(col="black") +
  geom_errorbar(aes(ymin=kmean_protect-kse_protect, ymax=kmean_protect+kse_protect))+
  geom_errorbarh(aes(xmin=rmean_protect-rse_protect, xmax=rmean_protect+rse_protect)) + 
  geom_label(aes(x=rmean_path, y=kmean_path, label=path))  +
  geom_point(aes(x=rmean_path, y=kmean_path))  +
  geom_point(data=merged_isointerest, aes(col=protect), cex=3) +
  ylab("Carrying capacity (k)") + xlab("Growth rate (r)")+
  labs(col="Protective Strain\nof interest")
gg_ecology_protect
ggsave("03_EDA/gg_ecology_protect.png", gg_ecology_protect, width=6, height=4)

######## Biofilm?
# merged$medianOD_path
merged %>%
  filter(path%in%c("ATR2A1","ATY2A2")) %>%
  ggplot(aes(x=rmean_protect, y=medianOD_protect)) + 
  geom_point(col="black") +
  # geom_errorbar(aes(ymin=kmean_protect-kse_protect, ymax=kmean_protect+kse_protect))+
  # geom_errorbarh(aes(xmin=rmean_protect-rse_protect, xmax=rmean_protect+rse_protect)) + 
  geom_label(aes(x=rmean_path, y=medianOD_path, label=path))  +
  geom_point(aes(x=rmean_path, y=medianOD_path))  +
  geom_point(data=merged_isointerest, aes(col=protect), cex=3) +
  ylab("Biofilm thickness (CV)") + xlab("Growth rate (r)")+
  labs(col="Protective Strain\nof interest")
# gg_ecology_protect
# ggsave("03_EDA/gg_ecology_protect.png", gg_ecology_protect, width=6, height=4)

gg_effect_of_biofilm <- merged %>%
  filter(!is.na(medianOD_protect)) %>%
  ggplot(aes(x=medianOD_protect, y=medianGA, col=path)) + 
  geom_point()+
  geom_smooth(method="lm")+
  ylab("Median green leaf area") + xlab("Median biofilm thickness\n(CV intensity)") +
  labs(col="Treatment") +
  scale_color_manual(values=c("red","orange","blue"))
gg_effect_of_biofilm
ggsave("03_EDA/gg_effect_of_biofilm.png", gg_effect_of_biofilm, width=6, height=4)

###### concentration gradient #####
# Grad_2023-02-01 is gradient of pathogen test
# Gt is separate

plant_adj %>%
  filter(experiment=="Grad_2023-02-01") %>%
  ggplot() +
  geom_boxplot(aes(x=factor(od_path), y=greenarea, col=path))+
  geom_point(aes(x=factor(od_path), y=greenarea, col=path, group=path), position = position_dodge(width=0.75))
  
plant_gt %>%
  mutate(Additive=gsub("NA-","",Treatment)) %>%
  mutate(Additive = factor(Additive, levels=c("MOCK","LB","glycerol","glycerolLB"))) %>%
  mutate(dilution = (1/3)^(col-1), dilution_factor = ifelse(col==1, 1, col-1)) %>%
  # select(col,dilution) %>%View()
  mutate(concentration = ifelse(col==1,0,dilution)) %>%
  # filter(col!=1) %>%
  ggplot(aes(x=dilution_factor, y=greenarea, col=Additive)) +
  geom_point() +
  geom_smooth(se = FALSE)
  # geom_point(aes(x=factor(od_path), y=greenarea, col=path, group=path), position = position_dodge(width=0.75))

gg_gteffects <- plant_gt %>%
  mutate(Additive=gsub("NA-","",Treatment)) %>%
  mutate(Additive = factor(Additive, levels=c("MOCK","LB","glycerol","glycerolLB"))) %>%
  mutate(dilution = (1/3)^(col-1), dilution_factor = ifelse(col==1, 1, col-1)) %>%
  # select(col,dilution) %>%View()
  mutate(concentration = ifelse(col==1,0,dilution)) %>%
  # filter(col!=1) %>%
  ggplot(aes(x=factor(dilution_factor), y=greenarea, col=Additive)) +
  geom_boxplot() +
geom_point(aes(x=factor(dilution_factor), y=greenarea, col=Additive, group=Additive), position = position_dodge(width=0.75)) +
  xlab("Dilution factor (1/3)^n") +
  ylab("Green leaf area (sqrt pixels)")

ggsave("03_EDA/gg_gt_effects.png", gg_gteffects, width=6, height=4)

# Check glycerol??