library(tidyverse)
library(lubridate)
library(questionr)

## Read in SARS2 ORF8 dataset
ko <- read_tsv("washington.orf8ko.tsv")
metadata <- read_tsv("data/2022-07-26.gisaid.washington.metadata.tsv")
wdrs <- read_csv("../data/wdrs_metadata_9.2.22.csv")


## Merge datasets
# First, let's see what would not join:
not_in_ko <- metadata %>%
  anti_join(ko, by="strain")
# For whatever reason, there are 17 samples not in ko file, which means they 
# were not in the alignment that went into the script. Probably got filtered out 
# in align step.
ko %>%
  anti_join(metadata,by="strain")
# But all the other samples in ko are in metadata...

not_in_wdrs <- metadata %>%
  anti_join(wdrs, by = c("strain" = "gisaid_id"))
not_in_wdrs %>%
  select(submitting_lab) %>%
  distinct() %>%
  print(n=100)
wa_not_in_wdrs <- not_in_wdrs %>%
  filter(str_detect(submitting_lab, "Washington|UW|Seattle|Fred Hutch|Atlas|Altius"))
# Hmm, there are 11,442 samples in metadata not in WDRS and definitely from 
# local submitters... This is more than the 5,000 or so Lauren mentioned. I'll 
# have to ask about it.
not_in_metadata_in_wdrs <- wdrs %>%
  anti_join(metadata, by = c("gisaid_id" = "strain"))
# There's 7564 samples in wdrs, not in metadata, I hope they are newer samples 
# added recently, but I have no way of checking this in the absence of metadata.
# Let's actually do the join
df <- metadata %>%
  left_join(ko, by = "strain") %>%
  inner_join(wdrs, by = c("strain" = "gisaid_id"))

## Call ORF 8 KOs
df <- df %>%
  mutate(orf8ko = ifelse(gap > 6 | protein_length <= 100 | Ns > 100, "yes", "no"))

## Plot distribution of potential predictors by orf8ko
# Sex at birth
df %>%
  ggplot(aes(orf8ko, color=sex_at_birth,fill=sex_at_birth)) +
  geom_bar(position="fill")
# Pretty even sex split in yes & no

df %>%
  ggplot(aes(orf8ko, color=sex_at_birth,fill=sex_at_birth)) +
  geom_bar(stat="count")
# And there are very few NA's so should be fiiine

# Age
df %>%
  ggplot(aes(orf8ko, color=age_group,fill=age_group)) +
  geom_bar(position="fill")
# There are younger people on average with the ko, than without. Not sure how much of that is due to variant/timing bias. But interesting...
# Age without alpha
df %>%
  filter(Nextstrain_clade != "20I (Alpha, V1)") %>%
  ggplot(aes(orf8ko, color=age_group,fill=age_group)) +
  geom_bar(position="fill")
## Age distributions essentially identical without alpha.

# Variant
df %>%
  ggplot(aes(orf8ko, color=Nextstrain_clade,fill=Nextstrain_clade)) +
  geom_bar(position="fill")
# Alpha is confusing this, since all of alpha have the orf8 ko. I'll compare without alpha.
df %>%
  filter(Nextstrain_clade != "20I (Alpha, V1)") %>%
  ggplot(aes(orf8ko, color=Nextstrain_clade,fill=Nextstrain_clade)) +
  geom_bar(position="fill")
# without Alpha we do not have the same distribution across variants in orf8 KOs.
# This is a good reason to stick with the original plan of matched controls.

## Vaccines
df <- df %>%
  mutate(IIS_VACCINE_INFORMATION_AVAILABLE_DATE_1 = as_date(ifelse(difftime(date,IIS_VACCINE_INFORMATION_AVAILABLE_DATE_1, units="days") < 14, NA, IIS_VACCINE_INFORMATION_AVAILABLE_DATE_1))) %>%
  mutate(IIS_VACCINE_INFORMATION_AVAILABLE_DATE_2 = as_date(ifelse(difftime(date,IIS_VACCINE_INFORMATION_AVAILABLE_DATE_2, units="days") < 14, NA, IIS_VACCINE_INFORMATION_AVAILABLE_DATE_2))) %>%
  mutate(IIS_VACCINE_INFORMATION_AVAILABLE_DATE_3 = as_date(ifelse(difftime(date,IIS_VACCINE_INFORMATION_AVAILABLE_DATE_3, units="days") < 14, NA, IIS_VACCINE_INFORMATION_AVAILABLE_DATE_3))) %>%
  mutate(IIS_VACCINE_INFORMATION_AVAILABLE_DATE_4 = as_date(ifelse(difftime(date,IIS_VACCINE_INFORMATION_AVAILABLE_DATE_4, units="days") < 14, NA, IIS_VACCINE_INFORMATION_AVAILABLE_DATE_4))) %>%
  mutate(last_vaccination = pmax(IIS_VACCINE_INFORMATION_AVAILABLE_DATE_1,IIS_VACCINE_INFORMATION_AVAILABLE_DATE_2, IIS_VACCINE_INFORMATION_AVAILABLE_DATE_3, IIS_VACCINE_INFORMATION_AVAILABLE_DATE_4, na.rm=TRUE)) %>%
  mutate(days_since_vaccination = difftime(date, last_vaccination, units="days")) %>%
  mutate(vaccinated = ifelse(is.na(days_since_vaccination), "no", "yes"))
# Vaccine is eliminated if received less than 14 days ago. After elimination, 
# calculates day since last vaccine was given.
df %>%
  filter(Nextstrain_clade != "20I (Alpha, V1)") %>%
  filter(!is.na(orf8ko)) %>%
  ggplot(aes(days_since_vaccination,color=orf8ko,fill=orf8ko,)) +
  geom_histogram(aes(y=0.5*..density..),alpha=0.5,position='identity')
# Distributions seem to be very similar for days_since_vaccination by orf8ko

## Hospital
df$hosp = as_factor(df$hosp)
df$hosp = fct_relevel(df$hosp, "No", "Yes")
df %>%
  filter(Nextstrain_clade != "20I (Alpha, V1)") %>%
  filter(!is.na(orf8ko)) %>%
  ggplot(aes(orf8ko, color=hosp,fill=hosp)) +
  geom_bar(position="fill")
# Doesn't look to be a lot of difference
# Which clades have more than 30 orf8kos?
clades_30 <- df %>%
  group_by(Nextstrain_clade, orf8ko) %>%
  summarise(number = n()) %>%
  arrange(Nextstrain_clade) %>%
  filter(number > 30) %>%
  group_by(Nextstrain_clade) %>%
  summarise(rows = n()) %>%
  filter(rows > 1) %>%
  pull(Nextstrain_clade)

df %>%
  filter(Nextstrain_clade %in% clades_30) %>%
  filter(!is.na(orf8ko)) %>%
  filter(!is.na(hosp)) %>%
  ggplot(aes(hosp, color=orf8ko,fill=orf8ko)) +
  geom_bar(position="fill") +
  facet_wrap(~Nextstrain_clade)

df %>%
  filter(Nextstrain_clade %in% clades_30) %>%
  filter(!is.na(hosp)) %>%
  filter(!is.na(orf8ko)) %>%
  ggplot(aes(orf8ko, color=hosp,fill=hosp)) +
  geom_bar(position="fill") +
  facet_wrap(~Nextstrain_clade)

df %>%
  filter(Nextstrain_clade %in% clades_30) %>%
  filter(!is.na(hosp)) %>%
  filter(!is.na(orf8ko)) %>%
  filter(hosp != "Unknown") %>%
  ggplot(aes(orf8ko, color=hosp,fill=hosp)) +
  geom_bar(position="fill") +
  facet_wrap(~Nextstrain_clade)
# NOT much difference

## Death
df <- df %>%
  mutate(died = ifelse(is.na(death_date), "no", "yes"))

df %>%
  filter(Nextstrain_clade != "20I (Alpha, V1)") %>%
  filter(!is.na(orf8ko)) %>%
  ggplot(aes(orf8ko, color=died,fill=died)) +
  geom_bar(position="fill")

df %>%
  filter(Nextstrain_clade %in% clades_30) %>%
  filter(!is.na(orf8ko)) %>%
  ggplot(aes(orf8ko, color=died,fill=died)) +
  geom_bar(position="fill") +
  facet_wrap(~Nextstrain_clade)
# Doesn't look to be a lot of difference


## Logistic regression
# Hospitalization
df$age_group = fct_relevel(df$age_group, "0-4", "5-17", "18-44", "45-64", "65-79", "80+", "Unknown")
df_hosp <- df %>%
  filter(hosp != "Unknown") %>%
  filter(age_group != "Unknown") %>%
  filter(!is.na(age_group)) %>%
  filter(!is.na(orf8ko))
df_hosp$age_group = as.numeric(df_hosp$age_group)



clades_30_hosp <- df_hosp %>%
  group_by(Nextstrain_clade, orf8ko) %>%
  summarise(number = n()) %>%
  arrange(Nextstrain_clade) %>%
  filter(number > 30) %>%
  group_by(Nextstrain_clade) %>%
  summarise(rows = n()) %>%
  filter(rows > 1) %>%
  pull(Nextstrain_clade)


reg1 <- glm(hosp ~ orf8ko + age_group + sex_at_birth + vaccinated, family = "binomial", data = df_hosp)
summary(reg1)

questionr::odds.ratio(reg1, level=(1-(0.05)))


for (variant in clades_30_hosp) {
  reg <- glm(hosp ~ orf8ko + age_group + sex_at_birth + vaccinated, family = "binomial", data = df_hosp[df_hosp$Nextstrain_clade == variant,])
  summary(reg)
  or <- as.tibble(questionr::odds.ratio(reg, level=1-0.05),rownames="variable")
  print(variant)
  print(or)
}

# Death
df_death <- df %>%
  filter(age_group != "Unknown") %>%
  filter(!is.na(age_group)) %>%
  filter(!is.na(orf8ko)) %>%
  filter(!is.na(died))
df_death$age_group = as.numeric(df_death$age_group)
df_death$died = fct_relevel(df_death$died, "no", "yes")



clades_30_death <- df_death %>%
  group_by(Nextstrain_clade, orf8ko) %>%
  summarise(number = n()) %>%
  arrange(Nextstrain_clade) %>%
  filter(number > 30) %>%
  group_by(Nextstrain_clade) %>%
  summarise(rows = n()) %>%
  filter(rows > 1) %>%
  pull(Nextstrain_clade)


reg2 <- glm(died ~ orf8ko + age_group + sex_at_birth + vaccinated, family = "binomial", data = df_death)
summary(reg2)

questionr::odds.ratio(reg2, level=(1-(0.05)))


for (variant in clades_30_death) {
  reg <- glm(died ~ orf8ko + age_group + sex_at_birth + vaccinated, family = "binomial", data = df_death[df_death$Nextstrain_clade == variant,])
  summary(reg)
  or <- as.tibble(questionr::odds.ratio(reg, level=1-0.05),rownames="variable")
  print(variant)
  print(or)
}

