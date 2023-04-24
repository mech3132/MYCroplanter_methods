#!bin/bash
library(tidyverse)
library(brms)

dir.create("2022-06-17_statistics")

### load ####
dat <- read.delim("2022-06-14-analyze_data_so_far/dat_clean.txt") 
# Make Mock the plant value
allIsos <- unique(dat$isolate)
allIsos <- allIsos[allIsos!="MOCK"]

dat_forbrm_singles <- dat %>%
  mutate(Treatment = as.numeric(Type!="MOCK")) %>%
  mutate(isolate = factor(isolate, levels=c("MOCK", allIsos))) %>%
  select(Weight, Treatment, isolate, plate, experiment) %>%
  filter(experiment != "2022-06-07_rhizo_challenge") %>%
  drop_na() 
if (file.exists("2022-06-17_statistics/brm_iso_singles.RData")) {
  load("2022-06-17_statistics/brm_iso_singles.RData")
} else {
  brm_iso_singles <- brm(Weight ~ isolate +(1|experiment/plate), 
                 data   = dat_forbrm_singles, 
                 control = list(adapt_delta=0.99),
                 warmup = 1000, 
                 iter   = 4000, 
                 chains = 4) 
  save(brm_iso_singles, file="2022-06-17_statistics/brm_iso_singles.RData")
}


iso_draws <- as.data.frame(as_draws_df(brm_iso_singles))
iso_draws_b <- iso_draws %>% select(starts_with("b_isolate")) %>%
  rename_all(list(~gsub("b_isolate","",.))) %>%
  rownames_to_column(var="iter") %>%
  pivot_longer(-iter, names_to = "isolate", values_to = "samp") %>%
  # mutate(isolate = gsub("+O","PO",gsub("[+]F","PF",gsub("P","+",gsub("M","-",isolate))))) %>%
  group_by(isolate) %>% mutate(ave=mean(samp), sig = (quantile(samp, 0.005)>0 | quantile(samp, 0.995)<0)) %>% ungroup() %>%
  arrange(ave) %>% mutate(isolate = factor(isolate, levels=unique(isolate))) %>%
  select(-ave) %>%
  # mutate(isolate = factor(isolate, levels=c("MOCK", allIsos)))  %>%
  mutate(variable = "b") 
save(iso_draws_b, file = "2022-06-17_statistics/iso_draws_b.RData")

iso_draws_mock <- iso_draws %>% select("b_Intercept") %>%
  rename(MOCK = b_Intercept) %>%
  rownames_to_column(var="iter") %>%
  pivot_longer(-iter, names_to = "isolate", values_to = "samp") %>%
  mutate(variable = "Intercept")
save(iso_draws_mock, file = "2022-06-17_statistics/iso_draws_mock.RData")

iso_draws_randeff <- iso_draws %>% select(starts_with("r_experiment")) %>%
  # rename_all(funs(gsub("r_experiment","",.))) %>%
  rownames_to_column(var="iter") %>%
  pivot_longer(-iter, names_to = "isolate", values_to = "samp") %>%
  rowwise() %>%
  mutate(variable = ifelse(length(grep("plate", isolate))>0, "PlateEffect", "ExperimentEffect")) %>%
  mutate(experiment = gsub("_iso.*[,]Intercept[]]","_iso",gsub("r_experiment:plate[[]","",gsub("r_experiment[[]","",isolate)))) %>%
  mutate(plate = gsub(",Intercept[]]","",gsub("r_experiment:plate[[]2022.*iso_","",gsub("r_experiment[[].*_iso","",isolate)))) %>%
  ungroup() 
iso_draws_exp <- iso_draws_randeff %>% filter(variable=="ExperimentEffect") %>% 
  rename(exp_samp = samp) %>%
  select(iter, exp_samp, experiment) 
iso_draws_plt <- iso_draws_randeff %>% filter(variable=="PlateEffect") %>% 
  rename(plt_samp = samp ) %>%
  select(iter, plt_samp, experiment, plate)
iso_draws_allrandeff <- full_join(iso_draws_plt,iso_draws_exp) %>%
  mutate(samp = plt_samp + exp_samp)
save(iso_draws_allrandeff, file = "2022-06-17_statistics/iso_draws_allrandeff.RData")


## Raw data setup
iso_raw_dat_singles_adj <- iso_draws_mock %>% select(iter, samp) %>%
  rename(mock_samp = samp) %>% full_join(iso_draws_allrandeff) %>%
  mutate(samp_total_mock = samp + mock_samp) %>%
  group_by(experiment, plate) %>% summarise(mock_mean = mean(samp_total_mock)) %>%
  ungroup() %>% right_join(dat_forbrm_singles) %>%
  mutate(diff_weight = Weight-mock_mean)%>%
  select(experiment, plate, isolate, diff_weight) %>%
  filter(isolate!="MOCK") 
save(iso_raw_dat_singles_adj, file = "2022-06-17_statistics/iso_raw_dat_singles_adj.RData")

### COMBINE
iso_draws_SINGLES <- iso_draws_b %>% left_join(dat %>% select(isolate, Type, plate, experiment) %>% distinct()) %>%
  # left_join(iso_raw_dat_singles_adj) %>%
  # filter(experiment!="2022-06-07_rhizo_challenge") %>%
  group_by(isolate) %>% mutate(ave=mean(samp)) %>% ungroup() %>%
  arrange(ave) %>% mutate(isolate = factor(isolate, levels=unique(isolate))) %>%
  select(-ave)
save(iso_draws_SINGLES, file = "2022-06-17_statistics/iso_draws_SINGLES.RData")

gg_singleeffects <- iso_draws_SINGLES %>%
  ggplot() +
  geom_violin(aes(x=isolate, y=samp))+
  geom_point(dat=iso_raw_dat_singles_adj, aes(x=isolate, y=diff_weight), col="green") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  geom_hline(aes(yintercept=0), col="red", lty=2)+
  xlab("Isolate") + ylab("Difference in weight (mg)\n from mock (10mM MgSO4)\n[Model Posterior distribution]")
gg_singleeffects

ggsave(filename = "2022-06-17_statistics/singleeffects.png"
       , gg_singleeffects, width=12, height=3.5)

gg_randomeffects <- iso_draws_allrandeff %>%
  mutate(plate=factor(plate, levels=paste0("p",seq(1,35)))) %>%
  ggplot() +
  geom_violin(aes(x=plate, y=samp))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  geom_hline(aes(yintercept=0), col="red", lty=2) +
  facet_grid(.~experiment, drop=TRUE, scales = "free") +
  xlab("Plate") + ylab("Deviation from total average, mg [Model Posterior distribution]")
gg_randomeffects
ggsave(filename = "2022-06-17_statistics/randomeffects.png"
       , gg_randomeffects, width=10, height=3.5)




## Challenge ######
dat_forbrm_challenge <- dat %>%
  mutate(Treatment = as.numeric(Type!="MOCK")) %>%
  mutate(isolate = factor(isolate, levels=c("MOCK", allIsos))) %>%
  filter(experiment=="2022-06-07_rhizo_challenge") %>%
  select(Weight, Treatment, isolate, plate, experiment) %>%
  drop_na()
if (file.exists("2022-06-17_statistics/brm_iso_challenge.RData")) {
  load("2022-06-17_statistics/brm_iso_challenge.RData")
} else {
  brm_iso_challenge <- brm(Weight ~ isolate +(1|plate),
                 data   = dat_forbrm_challenge,
                 control = list(adapt_delta=0.99),
                 warmup = 1000,
                 iter   = 4000,
                 chains = 4)
  save(brm_iso_challenge, file="2022-06-17_statistics/brm_iso_challenge.RData")
}

iso_draws_chall <- as.data.frame(as_draws_df(brm_iso_challenge))
iso_draws_chall_b <- iso_draws_chall %>% select(starts_with("b_isolate")) %>%
  rename_all(list(~gsub("b_isolate","",.))) %>%
  rownames_to_column(var="iter") %>%
  pivot_longer(-iter, names_to = "isolate", values_to = "samp") %>%
  # rowwise() %>%
  # mutate(isolate = gsub("[+]O","PO",gsub("[+]F","PF",gsub("P","[+]",gsub("M","[-]",isolate))))) %>%
  # ungroup() %>%
  group_by(isolate) %>% mutate(ave=mean(samp)) %>% ungroup() %>%
  arrange(ave) %>% mutate(isolate = factor(isolate, levels=unique(isolate))) %>%
  select(-ave) %>%
  # mutate(isolate = factor(isolate, levels=c("MOCK", challIsos)))  %>%
  mutate(variable = "b")
save(iso_draws_chall_b, file = "2022-06-17_statistics/iso_draws_chall_b.RData")

iso_draws_chall_mock <- iso_draws_chall %>% select("b_Intercept") %>%
  rename(MOCK = b_Intercept) %>%
  rownames_to_column(var="iter") %>%
  pivot_longer(-iter, names_to = "isolate", values_to = "samp") %>%
  mutate(variable = "Intercept")
save(iso_draws_chall_mock, file = "2022-06-17_statistics/iso_draws_chall_mock.RData")

iso_draws_chall_randeff <- iso_draws_chall %>% select(starts_with("r_plate")) %>%
  # rename_chall(funs(gsub("r_experiment","",.))) %>%
  rownames_to_column(var="iter") %>%
  pivot_longer(-iter, names_to = "isolate", values_to = "samp") %>%
  rowwise() %>%
  mutate(plate=gsub(",Intercept[]]","",gsub("r_plate[[]","",isolate))) %>% ungroup() %>%
  mutate(variable = "PlateEffect")  %>%
  select(-isolate)
save(iso_draws_chall_randeff, file = "2022-06-17_statistics/iso_draws_chall_randeff.RData")

dattemp <- dat %>% select(isolate, Type, plate, experiment, pathogen, commensal, Treatment) %>% 
  mutate(Treatment = gsub("and","\nand\n", Treatment)) %>%
  filter(experiment=="2022-06-07_rhizo_challenge", isolate!="MOCK") %>%
  mutate(isolate = gsub("[+]","P",gsub("[-]","M",isolate))) %>%
  distinct()
iso_draws_chall_final <- iso_draws_chall_b %>% full_join(dattemp) 
save(iso_draws_chall_final, file = "2022-06-17_statistics/iso_draws_chall_final.RData")




gg_chall <- iso_draws_chall_final %>%
  # filter(experiment =="2022-06-07_rhizo_challenge") %>%
  # group_by(isolate) %>% mutate(ave=mean(samp)) %>% ungroup() %>%
  # arrange(ave) %>% mutate(isolate = factor(isolate, levels=unique(isolate))) %>%
  # select(-ave) %>%
  ggplot() +
  geom_violin(aes(x=Treatment, y=samp))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  geom_hline(aes(yintercept=0), col="red", lty=2)+
  facet_grid(commensal~pathogen, drop = TRUE, scales = "free_x")
gg_chall
ggsave(filename = "2022-06-17_statistics/iso_challenge.png", height=15, width=15
       ,gg_chall)

####### Filter for figure
unique(dattemp$commensal)
path_to_keep <- c("N2C3","PAO1","ATR2A1","ATY2A2")
comm_to_keep <- c("WCS365","CHAO","CH267","ATT3B5","ATL3A3","ATD3B2")
dattempAlt <- dattemp %>%
  filter(pathogen %in% path_to_keep, commensal %in% comm_to_keep) %>%
  left_join(iso_draws_chall_b) %>%
  mutate(pathogen=factor(pathogen, levels=path_to_keep), commensal = factor(commensal, levels=comm_to_keep)) %>%
  group_by(isolate, pathogen, commensal, Treatment) %>%
  summarise(Mean=mean(samp), lwrCI = quantile(samp, 0.025), upprCI = quantile(samp, 0.975)) %>%
  ungroup()

controldats <- dattempAlt %>%
  filter(Treatment %in% c("Pathogen only", "Commensal only")) 

gg_chall_filt <-dattempAlt %>%
  filter(!Treatment %in% c("Pathogen only", "Commensal only")) %>%
  ggplot() +
  # geom_boxplot(aes(x=Treatment, y=samp, fill=Treatment))+
  geom_point(aes(x="Combined", y=Mean, col=Treatment)) +
  geom_segment(aes(x="Combined", xend="Combined", y=lwrCI, yend=upprCI)) +
  geom_hline(data=controldats, aes(yintercept=Mean, group=Treatment, col=Treatment))+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  xlab("")+
  geom_hline(aes(yintercept=0), col="darkgreen", lty=2)+
  facet_grid(pathogen~commensal, drop = TRUE, switch = "x")+
  scale_color_manual(values=c("Commensal only"="blue","Pathogen only"="red", "Pathogen \nand\n Commensal"="grey"))+
  ylab("Difference in weight (mg)\nfrom mock (10mM MgSO4)\n[Posterior distribution]") 
gg_chall_filt
ggsave(filename = "2022-06-17_statistics/iso_challenge_filt.png", height=4, width=6
       ,gg_chall_filt)

# 
# ?as_draws
# as_draws(fixef(brm_iso))
# brm_iso$fit
