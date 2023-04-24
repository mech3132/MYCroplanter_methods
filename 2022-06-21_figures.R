#!bin/bash

library(tidyverse)
dir.create("2022-06-07_figures")
load("2022-06-17_statistics/iso_draws_SINGLES.RData")
load("2022-06-17_statistics/iso_draws_allrandeff.RData")
load("2022-06-17_statistics/iso_draws_mock.RData")
load("2022-06-17_statistics/iso_raw_dat_singles_adj.RData")

gg_singleeffects <- iso_draws_SINGLES %>%
  ggplot() + 
  geom_violin(aes(x=isolate, y=samp))+
  geom_point(aes(x=isolate, y=diff_weight), col="green") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  geom_hline(aes(yintercept=0), col="red", lty=2)+
  xlab("Isolate") + ylab("Posterior distribution")
gg_singleeffects

ggsave(filename = "2022-06-17_statistics/singleeffects.png"
       , gg_singleeffects, width=10, height=3.5)

gg_randomeffects <- iso_draws_allrandeff %>%
  mutate(plate=factor(plate, levels=paste0("p",seq(1,35)))) %>%
  ggplot() + 
  geom_violin(aes(x=plate, y=samp))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  geom_hline(aes(yintercept=0), col="red", lty=2) +
  facet_grid(.~experiment, drop=TRUE, scales = "free") +
  xlab("Plate") + ylab("Posterior distribution")
gg_randomeffects
ggsave(filename = "2022-06-17_statistics/randomeffects.png"
       , gg_randomeffects, width=10, height=3.5)


