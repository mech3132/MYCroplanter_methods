library(tidyverse)

dir.create("01d_biofilm_process")

#### Load ######
dat <- read.csv("2023-03-24_growthcurves/dat/cv_2023-03-28_formatted.csv")
meta <- read.csv("2023-03-24_growthcurves/meta/meta_biofilm_2023-03-24.csv")

allCV <- full_join(dat, meta) %>%
  filter(!is.na(isolate), is.na(remove))

mockOD <- allCV %>%
  filter(isolate=="MOCK") %>%
  group_by(experiment) %>%
  summarise(aveODMock = mean(OD))


# Manually remove some points because they are WAY outliers

allCV_edit <- allCV %>%
  left_join(mockOD) %>%
  mutate(OD_adj = OD - aveODMock) %>%
  filter(!(isolate=="ATL4A5" & OD_adj>0.5), !(isolate=="ATR3B4" & OD_adj>0.5)) 

allCV_summary <- allCV_edit%>%
  group_by(isolate) %>%
  summarise(medianOD = median(OD_adj, na.rm=TRUE), se = sd(OD_adj, na.rm=TRUE)/n())


gg_biofilm <- allCV_summary %>%
  arrange(medianOD) %>% mutate(isolate = factor(isolate, levels=c("MOCK",unique(isolate)[unique(isolate)!="MOCK"]))) %>%
  ggplot() +
  geom_point(aes(x=isolate, y=medianOD)) +
  geom_errorbar(aes(x=isolate, ymin=medianOD-se, ymax=medianOD+se)) +
  geom_point(data=allCV_edit, mapping=aes(x=isolate, y=OD_adj), col="yellow") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ylab("Median biofilm thickness\n(Crystal violet OD @ 590nm)")
gg_biofilm

ggsave(filename="01d_biofilm_process/biofilm_plot.png",
       gg_biofilm, width=6, height=3)


## Write out tables
write.table(allCV_edit, file="01d_biofilm_process/cv_data_full.txt", row.names = FALSE, quote = FALSE, sep="\t")
write.table(allCV_summary, file="01d_biofilm_process/cv_data_summary.txt", row.names = FALSE, quote = FALSE, sep="\t")

