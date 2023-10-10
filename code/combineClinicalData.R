library(tidyverse)
library(lubridate)

setwd("~/Work/projects/covid/long-deletions")
## Read in SARS2 ORF8 dataset
ko <- read_tsv("wa_results/gisaid.washington_ko_meta.tsv")
clustersAlpha <- read_tsv("nextstrain_build/results/Alpha/clusters/clusters_ORF8.tsv")
clustersDelta <- read_tsv("nextstrain_build/results/Delta/clusters/clusters_ORF8.tsv")
clustersOther <- read_tsv("nextstrain_build/results/WA_other/clusters/clusters_ORF8.tsv")
wdrs <- read_csv("data/wdrs_metadata_9.2.22.csv")

## Generate cluster size
Alpha_size <- clustersAlpha %>% group_by(cluster) %>%
  add_count() %>%
  ungroup %>%
  rename(clusterSize = n) %>%
  dplyr::select(strain,clusterSize)

Delta_size <- clustersDelta %>% group_by(cluster) %>%
  add_count() %>%
  ungroup %>%
  rename(clusterSize = n) %>%
  dplyr::select(strain,clusterSize)

Other_size <- clustersOther %>% group_by(cluster) %>%
  add_count() %>%
  ungroup %>%
  rename(clusterSize = n) %>%
  dplyr::select(strain,clusterSize)

size <- Alpha_size %>%
  bind_rows(Delta_size) %>%
  bind_rows(Other_size)


not_in_wdrs <- ko %>%
  anti_join(wdrs, by = c("strain" = "gisaid_id"))
not_in_wdrs %>%
  select(submitting_lab) %>%
  distinct() %>%
  print(n=100)
wa_not_in_wdrs <- not_in_wdrs %>%
  filter(str_detect(submitting_lab, "Washington|UW|Seattle|Fred Hutch|Atlas|Altius"))
# There are ~30,000 WA samples not in WDRS

wa_not_in_wdrs %>%
  ggplot(aes(x=date)) +
  geom_histogram()
## Many of these are after mid-2022 when the data was pulled. But there are ones
## missing from earlier. Oh well.



not_in_koin_wdrs <- wdrs %>%
  anti_join(ko, by = c("gisaid_id" = "strain"))
# There's 2152 samples in wdrs, not in metadata/ko. Maybe these were culled
# because there coverage was too low..

## Combine WDRS, cluster size, and KO info
df <- ko %>%
  inner_join(wdrs, by = c("strain" = "gisaid_id")) %>%
  left_join(size, by = c("strain")) %>%
  mutate(epiweek = epiweek(date)) %>%
  mutate(year = year(date)) %>%
  mutate(month = month(date)) %>%
  unite(wkyr, year,epiweek, remove=FALSE) %>%
  unite(moyr,month,year,remove=FALSE)

df$age_group = fct_relevel(df$age_group, "0-4", "5-17", "18-44", "45-64", "65-79", "80+")

# Vaccine is eliminated if received less than 14 days ago. After elimination,
# calculates day since last vaccine was given.
df <- df %>%
  mutate(IIS_VACCINE_INFORMATION_AVAILABLE_DATE_1 = as_date(ifelse(difftime(date,IIS_VACCINE_INFORMATION_AVAILABLE_DATE_1, units="days") < 14, NA, IIS_VACCINE_INFORMATION_AVAILABLE_DATE_1))) %>%
  mutate(IIS_VACCINE_INFORMATION_AVAILABLE_DATE_2 = as_date(ifelse(difftime(date,IIS_VACCINE_INFORMATION_AVAILABLE_DATE_2, units="days") < 14, NA, IIS_VACCINE_INFORMATION_AVAILABLE_DATE_2))) %>%
  mutate(IIS_VACCINE_INFORMATION_AVAILABLE_DATE_3 = as_date(ifelse(difftime(date,IIS_VACCINE_INFORMATION_AVAILABLE_DATE_3, units="days") < 14, NA, IIS_VACCINE_INFORMATION_AVAILABLE_DATE_3))) %>%
  mutate(IIS_VACCINE_INFORMATION_AVAILABLE_DATE_4 = as_date(ifelse(difftime(date,IIS_VACCINE_INFORMATION_AVAILABLE_DATE_4, units="days") < 14, NA, IIS_VACCINE_INFORMATION_AVAILABLE_DATE_4))) %>%
  mutate(last_vaccination = pmax(IIS_VACCINE_INFORMATION_AVAILABLE_DATE_1,IIS_VACCINE_INFORMATION_AVAILABLE_DATE_2, IIS_VACCINE_INFORMATION_AVAILABLE_DATE_3, IIS_VACCINE_INFORMATION_AVAILABLE_DATE_4, na.rm=TRUE)) %>%
  mutate(days_since_vaccination = difftime(date, last_vaccination, units="days")) %>%
  mutate(vaccinated = ifelse(is.na(days_since_vaccination), "no", "yes"))

df <- df %>%
  mutate(died = ifelse(is.na(death_date), "No", "Yes"))

df %>%
  select(hosp,age_group,ORF8_ko,sex_at_birth,vaccinated,days_since_vaccination, Nextstrain_clade,wkyr,clusterSize,died,moyr,coverage) %>%
  write_tsv('wa_results/clinical_combined.tsv')
