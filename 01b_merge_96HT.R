#!/bin/bash Rscript
library(tidyverse)
# library(lme4)
# library(lmerTest)
# library(brms)
dir.create("01b_merge_96HT")
### Merging two polymyxas ####

dat1 <- read.delim("2022-12-06_athal_96_polymyxan2c3/data/processed_image/pixeldata_green_ratio_mask.txt")
dat21 <- read.delim("2022-12-10_ahtal_96_polymyxan2c3/data/processed_athal_scan_001/pixeldata_green_ratio_mask.txt")
dat22 <- read.delim("2022-12-10_ahtal_96_polymyxan2c3/data/processed_athal_scan_002/pixeldata_green_ratio_mask.txt")
dat23 <- read.delim("2022-12-10_ahtal_96_polymyxan2c3/data/processed_athal_scan_003/pixeldata_green_ratio_mask.txt")
dat24 <- read.delim("2022-12-10_ahtal_96_polymyxan2c3/data/processed_athal_scan_004/pixeldata_green_ratio_mask.txt")
dat31 <- read.delim("2023-01-23_sar_gradientMESload/data/processed_scans/2023-01-31_scanned_mycgradientpath038/pixeldata_green_ratio_mask.txt")
dat32 <- read.delim("2023-01-23_sar_gradientMESload/data/processed_scans/2023-01-31_scanned_bottomsar002/pixeldata_green_ratio_mask.txt")
dat41 <- read.delim("2023-02-01_gradient_polymyxan2c3/data/processed_scans/polymyxa_plate002/pixeldata_green_ratio_mask.txt")
dat42 <- read.delim("2023-02-01_gradient_polymyxan2c3/data/processed_scans/N2C3_plate001/pixeldata_green_ratio_mask.txt")
dat51 <- read.delim("2023-03-08_polymyxa_screen/data/processed_scans/2023-03-15_scan001/pixeldata_green_ratio_mask.txt")
dat52 <- read.delim("2023-03-08_polymyxa_screen/data/processed_scans/2023-03-15_scan002/pixeldata_green_ratio_mask.txt")
dat61 <- read.delim("2023-03-15_polymyxa_screen_and_GT/data/processed_scans/athal_screen_polymyxa001/pixeldata_green_ratio_mask.txt")
dat62 <- read.delim("2023-03-15_polymyxa_screen_and_GT/data/processed_scans/athal_screen_polymyxa002/pixeldata_green_ratio_mask.txt")
dat63 <- read.delim("2023-03-15_polymyxa_screen_and_GT/data/processed_scans/athal_screen_polymyxa_gt003/pixeldata_green_ratio_mask.txt")
dat71 <- read.delim("2023-03-21_polymyxa_screen/data/processed_scans/athal_polymyxa_001/pixeldata_green_ratio_mask.txt")
dat72 <- read.delim("2023-03-21_polymyxa_screen/data/processed_scans/athal_polymyxa_002/pixeldata_green_ratio_mask.txt")
dat73 <- read.delim("2023-03-21_polymyxa_screen/data/processed_scans/athal_polymyxa_003/pixeldata_green_ratio_mask.txt")
dat74 <- read.delim("2023-03-21_polymyxa_screen/data/processed_scans/athal_polymyxa_004/pixeldata_green_ratio_mask.txt")


meta1 <- read.csv("2022-12-06_athal_96_polymyxan2c3/meta/meta_2022-12-06.csv")
meta2 <- read.csv("2022-12-10_ahtal_96_polymyxan2c3/meta/meta_2022-12-10.csv")
meta3 <- read.csv("2023-01-23_sar_gradientMESload/meta/meta_2023-01-23.csv")
meta4 <- read.csv("2023-02-01_gradient_polymyxan2c3/meta/meta_2023-02-01.csv")
meta5 <- read.csv("2023-03-08_polymyxa_screen/meta/meta_2023-03-08.csv")
meta6 <- read.csv("2023-03-15_polymyxa_screen_and_GT/meta/meta_2023-03-15.csv")
meta7 <- read.csv("2023-03-21_polymyxa_screen/meta/meta_2023-03-21.csv")

manremove1 <- read.delim("2022-12-06_athal_96_polymyxan2c3/data/manual_remove.txt", header=FALSE) %>%
  rename(subimageName = V1) %>% mutate(remove=TRUE)
manremove2 <- read.delim("2022-12-10_ahtal_96_polymyxan2c3/data/manual_remove.txt", header=FALSE) %>%
  rename(subimageName = V1, plate = V2) %>% mutate(remove=TRUE)
manremove3 <- read.csv("2023-01-23_sar_gradientMESload/data/manual_remove.txt", header=FALSE) %>%
  rename(subimageName = V1, plate = V2) %>% mutate(remove=TRUE)
manremove4 <- read.csv("2023-02-01_gradient_polymyxan2c3/data/manual_remove.txt", header=FALSE) %>%
  rename(subimageName = V1, plate = V2) %>% mutate(remove=TRUE)
manremove5 <- read.csv("2023-03-08_polymyxa_screen/data/manual_remove.txt", header=FALSE) %>%
  rename(subimageName = V1, plate = V2) %>% mutate(remove=TRUE)
manremove6 <- read.csv("2023-03-15_polymyxa_screen_and_GT/data/manual_remove.txt", header=FALSE) %>%
  rename(subimageName = V1, plate = V2) %>% mutate(remove=TRUE)
manremove7 <- read.csv("2023-03-21_polymyxa_screen/data/manual_remove.txt", header=FALSE) %>%
  rename(subimageName = V1, plate = V2) %>% mutate(remove=TRUE)

taxa <- read.csv("00_SANGERSEQUENCES/03_adding_taxonomy/sequence_id_taxonomy_MYConly.csv")

#### Merge ########
all1 <- dat1 %>% 
  left_join(manremove1) %>%
  filter(is.na(remove)) %>%
  mutate(experiment='Ht_2022-12-06', plate=1) %>%
  separate(subimageName, into=c("row","col"), sep = 1) %>%
  left_join(meta1 %>% mutate(col=as.character(col))) 

all2 <- dat21 %>% 
  mutate( plate=1) %>%
  full_join(dat22 %>% mutate(plate=2)) %>%
  full_join(dat23 %>% mutate(plate=3)) %>%
  full_join(dat24 %>% mutate(plate=4)) %>%
  left_join(manremove2) %>%
  filter(is.na(remove)) %>%
  mutate(experiment = 'Ht_2022-12-10') %>%
  separate(subimageName, into=c("row","col"), sep = 1) %>% 
  left_join(meta2 %>% mutate(col=as.character(col))) %>%
  rename(ratio_PR_TO_PA = ratio)

all3 <- dat31 %>% 
  mutate( plate=1) %>%
  full_join(dat32 %>% mutate(plate=2)) %>%
  # full_join(dat23 %>% mutate(plate=3)) %>%
  # full_join(dat24 %>% mutate(plate=4)) %>%
  left_join(manremove3) %>%
  filter(is.na(remove)) %>%
  mutate(experiment = 'Grad_2023-01-23') %>%
  separate(subimageName, into=c("row","col"), sep = 1) %>% 
  left_join(meta3 %>% mutate(col=as.character(col))) %>%
  rename(ratio_PR_TO_PA = ratio)

all4 <- dat41 %>% 
  mutate(plate=1) %>%
  full_join(dat42 %>% mutate(plate=2)) %>%
  left_join(manremove4) %>%
  filter(is.na(remove)) %>%
  mutate(experiment = 'Grad_2023-02-01') %>%
  separate(subimageName, into=c("row","col"), sep = 1) %>% 
  left_join(meta4 %>% mutate(col=as.character(col))) 


all5 <- dat51 %>% 
  mutate(plate=1) %>%
  full_join(dat52 %>% mutate(plate=2)) %>%
  left_join(manremove5) %>%
  filter(is.na(remove)) %>%
  mutate(experiment = 'Polymyxa_2023-03-08') %>%
  separate(subimageName, into=c("row","col"), sep = 1) %>% 
  left_join(meta5 %>% mutate(col=as.character(col))) %>%
  rename(ratio_PR_TO_PA = ratio)

all6 <- dat61 %>% 
  mutate(plate=1) %>%
  full_join(dat62 %>% mutate(plate=2)) %>%
  full_join(dat63 %>% mutate(plate=3)) %>%
  left_join(manremove6) %>%
  filter(is.na(remove)) %>%
  mutate(experiment = 'PolymyxaGT_2023-03-15') %>%
  separate(subimageName, into=c("row","col"), sep = 1) %>% 
  left_join(meta6 %>% mutate(col=as.character(col))) %>%
  rename(ratio_PR_TO_PA = ratio)

all7 <- dat71 %>% 
  mutate(plate=1) %>%
  full_join(dat72 %>% mutate(plate=2)) %>%
  full_join(dat73 %>% mutate(plate=3)) %>%
  full_join(dat73 %>% mutate(plate=4)) %>%
  left_join(manremove7) %>%
  filter(is.na(remove)) %>%
  mutate(experiment = 'Polymyxaother_2023-03-21') %>%
  separate(subimageName, into=c("row","col"), sep = 1) %>% 
  left_join(meta7 %>% mutate(col=as.character(col))) %>%
  rename(ratio_PR_TO_PA = ratio)

# all isolates
allIsos <- unique(c(unlist(all1[,c("path","protect")])
                  , unlist(all2[,c("path","protect")])
                  , unlist(all3[,c("path","protect")])
                  , unlist(all4[,c("path","protect")])
                  , unlist(all5[,c("path","protect")])
                  , unlist(all6[,c("path","protect")])
                  , unlist(all7[,c("path","protect")])))
allIsos <- allIsos[!allIsos%in%c("MOCK",NA, "glycerol")]

allProtect <- unique(c(unlist(all1[,c("protect")])
                    , unlist(all2[,c("protect")])
                    , unlist(all3[,c("protect")])
                    , unlist(all4[,c("protect")])
                    , unlist(all5[,c("protect")])
                    , unlist(all6[,c("protect")])
                    , unlist(all7[,c("protect")])))
allProtect <- allProtect[!allProtect%in%c("MOCK",NA)]

taxa_mapping <- taxa %>%
  filter(isolateID %in% allIsos) %>%
  select(isolateID, ScientificName) %>%
  filter(ScientificName != "Unassigned_sp.") %>% distinct()

mergedat <- full_join(all1, all2) %>%
  full_join(all3) %>%
  full_join(all4) %>%
  full_join(all5) %>%
  full_join(all6) %>%
  full_join(all7) %>%
  unite(protect, path, sep="-", col="Treatment", remove=FALSE) %>%
  mutate(path = factor(path, levels=c("MOCK","N2C3","ATR2A1","ATY2A2"))) %>%
  mutate(protect = factor(protect, levels=c("MOCK",allProtect))) %>%
  mutate(greenarea = sqrt(pixels)) %>%
  unite(experiment, plate, col="exp_plate", remove=FALSE) %>%
  mutate(plant_age_inoc = as.numeric(as.Date(date_inoc)-as.Date(date_germ))) %>%
  select(exp_plate, Treatment, protect, path, pixels, greenarea,row, col, experiment, plate, ratio_PR_TO_PA, od_path, od_protect,plant_age_inoc) %>%
  left_join(taxa_mapping %>% rename(protect=isolateID), multiple="all") %>%
  mutate(threshDA = greenarea>48) %>%
  mutate(clusterID = kmeans(greenarea, centers = 2)$cluster)%>%
  mutate(clusteraverage = mean(greenarea)) %>%
  mutate(clusterDA = ifelse(clusteraverage<50 | greenarea<45, 0, 1))

### Save file ####
write.table(mergedat, file="01b_merge_96HT/merged_plant_data.txt", row.names = FALSE, quote = FALSE, sep="\t")
write.table(taxa_mapping, file="01b_merge_96HT/taxa_plant_data.txt", row.names = FALSE, quote = FALSE, sep="\t")

# Plot results
## Bimodal
mergedat %>%
  ggplot() +
  geom_histogram(aes(x=greenarea, group=clusterDA, fill=clusterDA))+
  facet_wrap(.~experiment)

mergedat %>%
  ggplot() +
  geom_histogram(aes(x=greenarea, group=plant_age_inoc, fill=plant_age_inoc))+
  facet_wrap(.~experiment)

mergedat %>%
  filter(Treatment=="MOCK-MOCK") %>%
  # filter(experiment != "Ht_2022-12-10") %>%
  ggplot() +
  geom_point(aes(x=experiment, y=greenarea, col=as.factor(plant_age_inoc), pch=factor(plate)))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

mergedat %>%
  filter(experiment=="Grad_2023-02-01") %>%
  ggplot(aes(x=(od_path), y=greenarea)) + 
  geom_point()+
  geom_smooth(method="lm") +
  geom_hline(aes(yintercept=50), col="red") +
  facet_grid(.~od_protect) +
  xlab("OD of WCS365") + ylab("Green pixel area")

gg_ratios <- mergedat %>%
  filter(experiment=="Grad_2023-02-01") %>%
  ggplot(aes(x=(od_path), y=clusterDA)) +
  geom_jitter(height=0.1, width=0)+
  geom_smooth(method="glm",method.args = list(family = "binomial"))+
  geom_hline(aes(yintercept=0.5), col="red", lty=2)+
  facet_grid(.~od_protect) +
  xlab("OD of pathogen (N2C3)") 
gg_ratios

mergedat %>%
  filter(experiment=="Grad_2023-02-01") %>%
  ggplot(aes(x=(od_protect), y=clusterDA)) +
  geom_jitter(height=0.1, width=0)+
  geom_smooth(method="glm",method.args = list(family = "binomial"))+
  geom_hline(aes(yintercept=0.5), col="red")+
  facet_grid(.~od_path) +
  xlab("OD of protective (WCS365)") 

mergedat %>%
  filter(experiment=="Grad_2023-02-01") %>%
  ggplot(aes(x=(od_path), y=clusterDA, group=factor(log(od_protect+1)), col=log(od_protect+1))) +
  geom_jitter(height=0.1, width=0)+
  geom_smooth(method="glm",method.args = list(family = "binomial"), se = FALSE)+
  geom_hline(aes(yintercept=0.5), col="red")+
  # facet_grid(.~od_path) +
  xlab("OD of N2C3")
## Fit binomial model?
grad_filtset <- mergedat %>%
  filter(experiment=="Grad_2023-02-01") %>%
  mutate(ratioPRPA = od_protect/od_path)
ratio_glm <- glm(clusterDA ~ od_path*factor(od_protect), data=grad_filtset)
summary(ratio_glm)
ratio_glm2 <- glm(clusterDA ~ od_protect*factor(od_path), data=grad_filtset)
summary(ratio_glm2)
ratio_glm3 <- glm(clusterDA ~ od_protect*od_path, data=grad_filtset)
summary(ratio_glm3)

onlyGrad <- mergedat %>%
  filter(experiment=="Grad_2023-02-01")

##### model fit ######
hist(log(mergedat$pixels))
hist(log(mergedat$greenarea))
# Should log pixels



# Filter data to include only rhizosphere data
mergedat_filt <- mergedat %>%
  filter(protect %in% c(allIsos,"MOCK") | path %in% c(allIsos,"MOCK"))
# make a glycerol one
mergedat_gt <- mergedat %>%
  filter(experiment=="PolymyxaGT_2023-03-15", is.na(protect))
# I think I want to normalize pixel count by mock average
pix_mock <- mergedat_filt %>%
  filter(path=="MOCK"& protect=="MOCK") %>%
  group_by(experiment, plate) %>%
  summarise(meanga_mock = mean(greenarea),meanpix_mock = mean(pixels))
mergedat_filt <- mergedat_filt %>%
  filter(ratio_PR_TO_PA == "1 to 1") %>%
  left_join(pix_mock) %>%
  mutate(adjGA = (greenarea)/meanga_mock, adjPix = (pixels)/meanpix_mock) 

#### PLOTTING ########
# Just singles
# unique(mergedat_filt$experiment)




mergedat_filt %>%
  filter(path=="MOCK", !is.na(protect), !is.na(path)) %>%
  ggplot() + 
  geom_boxplot(aes(x=protect, y=adjGA)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
mergedat_filt %>%
  filter(path=="MOCK", !is.na(protect), !is.na(path)) %>%
  ggplot() + 
  geom_boxplot(aes(x=protect, y=adjPix, col=experiment)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

mergedat_filt %>%
  filter(protect=="MOCK", !is.na(protect), !is.na(path)) %>%
  ggplot(aes(x=path, y=adjGA)) + 
  geom_boxplot() +
  geom_jitter(width=0.25) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  ylab("Normalized pixels\n(:Mock)") +xlab("Pathogen")

mergedat_filt %>%
  filter(path%in% c("ATR2A1","MOCK"), !is.na(protect), !is.na(path)) %>%
  ggplot(aes(x=protect, y=adjGA, col=path)) + 
  geom_boxplot() +
  geom_jitter(width=0.25, height=0) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  ylab("Normalized pixels\n(:Mock)") 

mergedat_filt %>%
  filter(path%in% c("ATY2A2","MOCK"), !is.na(protect), !is.na(path)) %>%
  ggplot(aes(x=protect, y=adjGA, col=path)) + 
  geom_boxplot() +
  geom_jitter(width=0.25, height=0) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  ylab("Normalized pixels\n(:Mock)") 

taxa_mapping

# mergedat_filt %>% filter(adjPix==0)

lme_mdl <- lmer(log(adjPix) ~ path*protect + (1|experiment), data=mergedat_filt)
res_lme <- summary(lme_mdl) 

mdl_res_adj <- res_lme$coefficients %>%
  as.data.frame() %>%
  rownames_to_column(var="coef") %>%
  separate(coef, into=c("group1","group2"), sep=":", remove=FALSE) %>%
  rowwise() %>%
  mutate(path = ifelse( length(grep("path", coef)), group1, "MOCK")
         , protect = ifelse(is.na(group2) & length(grep("protect", coef)), group1, ifelse(is.na(group2),"MOCK",group2))) %>%
  mutate(path = gsub("path","",path), protect = gsub("protect","",protect)) %>%
  ungroup() %>%
  rename(se = `Std. Error`, p = `Pr(>|t|)`) %>%
  mutate(path=factor(path, levels=c("MOCK","N2C3","ATR2A1","ATY2A2"))
         , protect = factor(protect, levels=c("MOCK", allProtect)))

mdl_res_adj %>%
  ggplot(aes(x=protect)) +
  geom_bar(aes(y=Estimate, fill=path), stat="identity", position = position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=Estimate-se, ymax=Estimate+se,col=path), position = position_dodge(width=1))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

# if ( file.exists("01b_merge_96HT/bayes_mdl.RData")) {
#   load("01b_merge_96HT/bayes_mdl.RData")
# } else {
#   bayes_mdl <- brm(adjPix ~ path*protect + (1|experiment), data=mergedat_filt,family = lognormal)
#   save(bayes_mdl, file ="01b_merge_96HT/bayes_mdl.RData" )
# }


# #### Analysis #####
# 
# gg_byprotect <- mergedat %>%
#   filter(!is.na(path), !is.na(protect)) %>%
#   ggplot(aes(x=protect, y=sqrt(pixels))) + 
#   geom_boxplot() +
#   geom_jitter(cex=0.5, width=0.5, height=0) +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
#   facet_grid(experiment~path, scales="free", space = "free_x")
#   # xlab("Protective strain") +
#   # labs(title="Candidate protector used:") +
#   # ylab("Approximate leaf area\n(sqrt number green pixels)") +
#   # scale_color_manual(values=c("green","darkgreen"))
# gg_byprotect
# ############ MESSY 
# 
# ggsave("HT96_facet_by_protect.png"
#        , gg_byprotect
#        , width=7.5, height=3)
# 
# gg_bypathogen <- mergedat %>%
#   ggplot(aes(x=protect, y=sqrt(pixels), col=experiment)) + 
#   geom_boxplot() +
#   geom_jitter(cex=0.5, width=0.5, height=0) +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
#   facet_grid(.~path, scales = "free_x", drop=TRUE)+
#   xlab("Candidate protector strain") +
#   labs(title="Pathogen used to challenge:") +
#   ylab("Approximate leaf area\n(sqrt number green pixels)") +
#   scale_color_manual(values=c("green","darkgreen"))
# gg_bypathogen
# 
# ggsave("HT96_facet_by_pathogen.png"
#        , gg_bypathogen
#        , width=6, height=3)
# 
# ######## By exp plate
# 
# mergedat %>%
#   unite(experiment, plate, col="exp_plate", remove=FALSE) %>%
#   ggplot(aes(x=protect, y=sqrt(pixels), col=exp_plate)) + 
#   geom_boxplot() +
#   geom_jitter(cex=0.5, width=0.5, height=0) +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
#   facet_grid(exp_plate~path, scales = "free_x", drop=TRUE)+
#   xlab("Candidate protector strain") +
#   labs(title="Pathogen used to challenge:") +
#   ylab("Approximate leaf area\n(sqrt number green pixels)") 
# 
# mergedat %>%
#   unite(experiment, plate, col="exp_plate", remove=FALSE) %>%
#   ggplot(aes(x=path, y=sqrt(pixels), col=exp_plate)) + 
#   geom_boxplot() +
#   geom_jitter(cex=0.5, width=0.5, height=0) +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
#   facet_grid(exp_plate~protect, scales = "free_x", drop=TRUE)+
#   xlab("Pathogen strain") +
#   labs(title="Candidate protective used:") +
#   ylab("Approximate leaf area\n(sqrt number green pixels)") 
# 
# 
# 
# gg_summarypaeni <- mergedat %>%
#   unite(experiment, plate, col="exp_plate", remove=FALSE) %>%
#   filter(protect !="WCS365", path!="N2C3") %>%
#   filter(protect !="ATL4B3") %>%
#   ggplot(aes(x=path, y=sqrt(pixels), col=path)) + 
#   geom_boxplot(show.legend = FALSE) +
#   # geom_jitter(aes(col=exp_plate), cex=0.5, width=0.5, height=0) +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
#   facet_grid(.~protect, scales = "free_x", drop=TRUE)+
#   xlab("Pathogen strain") +
#   labs(title="Candidate protective used:") +
#   ylab("Approximate leaf area\n(sqrt number green pixels)") +
#   scale_color_manual(values=c("darkgreen","orange","darkorange"))
# gg_summarypaeni
# ggsave("HT96_paenionly.png"
#        ,gg_summarypaeni
#        , height=3, width=6)
# 
# 
# 
# 
# mergedat_filt <- mergedat %>%
#   filter(protect !="ATL4B3")
# 
# lm_basic <- lm(greenarea ~ path*protect, data = mergedat_filt)
# 
# lmer_basic <- lmer(greenarea ~ path*protect + (1|exp_plate) , data = mergedat_filt)
# # lmer_basic2 <- lmer(greenarea ~  (1+ path*protect|exp_plate) , data = mergedat_filt)
# 
# summary(lmer_basic)
# # summary(lmer_basic2)
# 
# CI_lmer <- confint(lmer_basic)
# CI_lmer %>%as.data.frame() %>% rownames_to_column(var="coef") %>%
#   ggplot() +
#   geom_segment(aes(x=coef, xend=coef, y=`2.5 %`, yend=`97.5 %`)) +
#   theme(axis.text.x = element_text(angle=90))
# 
# ########## Gradient dat ############
# all3 %>%
#   ggplot() + 
#   geom_point(aes(x=od_path, y=(pixels)))+
#   facet_wrap(.~path)
