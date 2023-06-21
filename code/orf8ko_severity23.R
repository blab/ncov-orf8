library(tidyverse)
library(lubridate)
library(questionr)

setwd("~/Work/projects/covid/long-deletions")
## Read in SARS2 ORF8 dataset
ko <- read_tsv("results/gisaid.washington.june20-july22.ORF8_ko.meta.tsv")
clusters <- read_tsv("clusters/detailedClusters/orf8ko_noAlpha_ko_meta_clusters.tsv")
wdrs <- read_csv("data/wdrs_metadata_9.2.22.csv")


## Merge datasets
# First, let's see what would not join:
not_in_ko <- metadata %>%
  anti_join(ko, by="strain")
ko %>%
  anti_join(metadata,by="strain")
# Everything is there

not_in_wdrs <- metadata %>%
  anti_join(wdrs, by = c("strain" = "gisaid_id"))
not_in_wdrs %>%
  select(submitting_lab) %>%
  distinct() %>%
  print(n=100)
wa_not_in_wdrs <- not_in_wdrs %>%
  filter(str_detect(submitting_lab, "Washington|UW|Seattle|Fred Hutch|Atlas|Altius"))
# Hmm, there are 9593 samples in metadata not in WDRS and definitely from 
# local submitters... This is more than the 5,000 or so Lauren mentioned. I'll 
# have to ask about it.
not_in_metadata_in_wdrs <- wdrs %>%
  anti_join(metadata, by = c("gisaid_id" = "strain"))
# There's 2142 samples in wdrs, not in metadata, I hope they are newer samples 
# added recently, but I have no way of checking this in the absence of metadata.
# Let's actually do the join

clustStatus = clusters %>%
  group_by(cluster) %>%
  count() %>% 
  rename(clusterSize=n) %>%
  right_join(clusters, by='cluster') %>%
  select(strain,cluster,clusterSize)

df <- ko %>%
  inner_join(wdrs, by = c("strain" = "gisaid_id")) %>%
  left_join(clustStatus) %>%
  mutate(koConfirm = ifelse(clusterSize>1|Nextstrain_clade=="20I (Alpha, V1)",'yes','no')) %>%
  mutate(koConfirm = ifelse(is.na(koConfirm),'no',koConfirm))

## Call ORF 8 KOs
df <- df %>%
  mutate(Nextstrain_clade = as_factor(Nextstrain_clade))

## Plot distribution of potential predictors by ORF8_ko
# Sex at birth
df %>%
  ggplot(aes(koConfirm, color=sex_at_birth,fill=sex_at_birth)) +
  geom_bar(position="fill")
# Pretty even sex split in yes & no

df %>%
  ggplot(aes(koConfirm, color=sex_at_birth,fill=sex_at_birth)) +
  geom_bar(stat="count")
# And there are very few NA's so should be fiiine

# Age
df %>%
  ggplot(aes(koConfirm, color=age_group,fill=age_group)) +
  geom_bar(position="fill")
# There are younger people on average with the ko, than without. Not sure how much of that is due to variant/timing bias. But interesting...
# Age without alpha
df %>%
  filter(Nextstrain_clade != "20I (Alpha, V1)") %>%
  ggplot(aes(ORF8_ko, color=age_group,fill=age_group)) +
  geom_bar(position="fill")
## Age distributions essentially identical without alpha.

# Variant
df %>%
  ggplot(aes(koConfirm, color=Nextstrain_clade,fill=Nextstrain_clade)) +
  geom_bar(position="fill")
# Alpha is confusing this, since all of alpha have the orf8 ko. I'll compare without alpha.
df %>%
  filter(Nextstrain_clade != "20I (Alpha, V1)") %>%
  ggplot(aes(ORF8_ko, color=Nextstrain_clade,fill=Nextstrain_clade)) +
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
  filter(!is.na(ORF8_ko)) %>%
  ggplot(aes(days_since_vaccination,color=ORF8_ko,fill=ORF8_ko,)) +
  geom_histogram(aes(y=0.5*..density..),alpha=0.5,position='identity')
# Distributions seem to be very similar for days_since_vaccination by ORF8_ko

## Hospital
df$hosp = as_factor(df$hosp)
df$hosp = fct_relevel(df$hosp, "No", "Yes","Unknown")
df %>%
  #filter(hosp!='Unknown') %>%
  filter(!is.na(ORF8_ko)) %>%
  filter(!is.na(hosp)) %>%
  ggplot(aes(koConfirm, color=hosp,fill=hosp)) +
  geom_bar(position="fill") +
  theme_minimal() +
  scale_x_discrete(labels=c('No','Yes')) +
  xlab('ORF8 KO') +
  ylab('Proportion') +
  scale_fill_manual(values=c('#457b9d','#14213d','#8d99ae'),name='Hospitalized?') +
  scale_color_manual(values=c('#457b9d','#14213d','#8d99ae'),name='Hospitalized?') 
  #scale_fill_manual(values=c('#457b9d','#14213d'),name='Hospitalized?') +
  #scale_color_manual(values=c('#457b9d','#14213d'),name='Hospitalized?')
ggsave('figs/hosp_bars.pdf',dpi=300,height=3,width=3)
# Doesn't look to be a lot of difference


## Death
df <- df %>%
  mutate(died = ifelse(is.na(death_date), "No", "Yes"))

df$death = as_factor(df$death)
df$death = fct_relevel(df$death, "No", "Yes")

df %>%
  filter(!is.na(ORF8_ko)) %>%
  filter(!is.na(died)) %>%
  ggplot(aes(koConfirm, color=died,fill=died)) +
  geom_bar(position="fill") +
  scale_x_discrete(labels=c("No","Yes")) +
  theme_minimal() +
  xlab('ORF8 KO') +
  ylab('Proportion') +
  scale_fill_manual(values=c('#457b9d','#14213d'),name='Death due to\nSARS-CoV-2?') +
  scale_color_manual(values=c('#457b9d','#14213d'),name='Death due to\nSARS-CoV-2?') 
#scale_fill_manual(values=c('#457b9d','#14213d'),name='Hospitalized?') +
#scale_color_manual(values=c('#457b9d','#14213d'),name='Hospitalized?')
ggsave('figs/death_bars.pdf',dpi=300,height=3,width=3)




## Logistic regression
# Hospitalization
df$age_group = fct_relevel(df$age_group, "0-4", "5-17", "18-44", "45-64", "65-79", "80+", "Unknown")
df_hosp <- df %>%
  filter(hosp != "Unknown") %>%
  filter(age_group != "Unknown") %>%
  filter(!is.na(age_group)) %>%
  filter(!is.na(ORF8_ko))
df_hosp$age_group = as.numeric(df_hosp$age_group)
df_hosp$hosp = as.factor(df_hosp$hosp)



clades_30_hosp <- df_hosp %>%
  group_by(Nextstrain_clade, ORF8_ko) %>%
  summarise(number = n()) %>%
  arrange(Nextstrain_clade) %>%
  filter(number > 30) %>%
  group_by(Nextstrain_clade) %>%
  summarise(rows = n()) %>%
  filter(rows > 1) %>%
  pull(Nextstrain_clade)


reg1 <- glm(hosp ~ koConfirm+ age_group + sex_at_birth + vaccinated, family = "binomial", data = df_hosp)
summary(reg1)

or_hosp = questionr::odds.ratio(reg1, level=(1-(0.05)))


for (variant in clades_30_hosp) {
  reg <- glm(hosp ~ ORF8_ko + age_group + sex_at_birth + vaccinated, family = "binomial", data = df_hosp[df_hosp$Nextstrain_clade == variant,])
  summary(reg)
  or <- as.tibble(questionr::odds.ratio(reg, level=1-0.05),rownames="variable")
  print(variant)
  print(or)
}

# Death
df_death <- df %>%
  filter(age_group != "Unknown") %>%
  filter(!is.na(age_group)) %>%
  filter(!is.na(ORF8_ko)) %>%
  filter(!is.na(died))
df_death$age_group = as.numeric(df_death$age_group)
df_death$died = fct_relevel(df_death$died, "no", "yes")



clades_30_death <- df_death %>%
  group_by(Nextstrain_clade, ORF8_ko) %>%
  summarise(number = n()) %>%
  arrange(Nextstrain_clade) %>%
  filter(number > 30) %>%
  group_by(Nextstrain_clade) %>%
  summarise(rows = n()) %>%
  filter(rows > 1) %>%
  pull(Nextstrain_clade)


reg2 <- glm(died ~ koConfirm + age_group + vaccinated, family = "binomial", data = df_death)
summary(reg2)

or_death = questionr::odds.ratio(reg2, level=(1-(0.05)))
or_death <- cbind(Variable = rownames(or_death), or_death)
or_hosp <- cbind(Variable = rownames(or_hosp), or_hosp)
or_hosp <- cbind(ko = c('no','yes','no','no','no','no'), or_hosp)

or_death <- cbind(ko = c('no','yes','no','no'), or_death)
colnames(or_hosp)[3] ="minCI"
colnames(or_hosp)[4] ="maxCI"
colnames(or_death)[3] ="minCI"
colnames(or_death)[4] ="maxCI"


or_hosp %>%
  filter(Variable!="(Intercept)") %>%
  ggplot(aes(y=Variable,x=OR)) +
  geom_vline(xintercept=1,linetype='dashed') +
  geom_linerange(aes(xmin=minCI,xmax=maxCI)) +
  geom_point(size=5,aes(color=ko,fill=ko))+
  scale_x_log10() +
  theme_minimal() +
  scale_y_discrete(labels=c('Age group', 'ORF8 KO','Sex: Male', 'Sex: Other', 'Vaccinated')) +
  ylab('') +
  xlab('Odds ratio for hospitalization') +
  scale_color_manual(values=c('grey','#ce4257')) +
  scale_fill_manual(values=c('grey','#ce4257')) +
  theme(legend.position = 'none')
ggsave('figs/or_hosp.pdf',dpi=300,height=3,width=4)

or_death %>%
  filter(Variable!="(Intercept)") %>%
  ggplot(aes(y=Variable,x=OR)) +
  
  geom_vline(xintercept=1,linetype='dashed') +
  geom_linerange(aes(xmin=minCI,xmax=maxCI)) +
  geom_point(size=5,aes(color=ko,fill=ko))+
  scale_x_log10() +
  theme_minimal() +
  scale_y_discrete(labels=c('Age group', 'ORF8 KO','Vaccinated')) +
  ylab('') +
  xlab('Odds ratio for death') +
  scale_color_manual(values=c('grey','#ce4257')) +
  scale_fill_manual(values=c('grey','#ce4257')) +
  theme(legend.position = 'none')
ggsave('figs/or_death.pdf',dpi=300,height=3,width=4)
                       






## Broken up by groups: Delta, Omicron, Gamma & 20A|20B|20C|20D|20G
group = list("Delta", "Omicron", "Gamma")


for (variant in group){
  cat(paste("##",variant,"\n","*Hospitalization:*\n"))
  freq_tab_hosp <- df_hosp %>%
    filter(variant %in% Nextstrain_clade) %>%
    group_by(ORF8_ko) %>%
    summarise(n_samples = n())
  print(kable(freq_tab_hosp, caption = paste("Samples with hospitalization data in ", variant, ".")))
  if (variant %in% clades_30_hosp){
    reg_hosp <- glm(hosp ~ ORF8_ko + age_group + sex_at_birth + vaccinated, family = "binomial", data = df_hosp[variant %in% df_hosp$Nextstrain_clade,])
    or_hosp <- as.tibble(questionr::odds.ratio(reg_hosp, level=(1-(0.05))), rownames = "Predictor")
    plot_hosp <- df_hosp %>%
      filter(variant %in% Nextstrain_clade) %>%
      ggplot(aes(ORF8_ko, color=hosp,fill=hosp)) +
      geom_bar(position="fill",width=0.7) +
      ylab("Proportion") +
      xlab("ORF8_ko") + 
      scale_fill_manual(values =c("#ADEFD1FF","#00203FFF")) +
      scale_color_manual(values =c("#ADEFD1FF","#00203FFF")) +
      theme_minimal() +
      theme(text = element_text(size = 10)) +
      ggtitle(label=variant)
    print(plot_hosp)
    print(kable(or_hosp[2:6,], caption = paste("Odds Ratio of Hospitalization for ",variant,"."),align="lllll"))
  }
  cat(paste("*Death:*\n"))
  freq_tab_death <- df_death %>%
    filter(variant %in% Nextstrain_clade) %>%
    group_by(ORF8_ko) %>%
    summarise(n_samples = n())
  print(kable(freq_tab_death, caption = paste("Samples with death data in ", variant, ".")))
  reg_death <- glm(died ~ ORF8_ko + age_group + sex_at_birth + vaccinated, family = "binomial", data = df_death[variant %in% df_hosp$Nextstrain_clade,])
  or_death <- as.tibble(questionr::odds.ratio(reg_death, level=(1-(0.05))), rownames = "Predictor")
  plot_death <- df_death %>%
    filter(variant %in% df_hosp$Nextstrain_clade) %>%
    ggplot(aes(ORF8_ko, color=died,fill=died)) +
    geom_bar(position="fill",width=0.7) +
    ylab("Proportion") +
    xlab("ORF8_ko") + 
    scale_fill_manual(values =c("#ADEFD1FF","#00203FFF")) +
    scale_color_manual(values =c("#ADEFD1FF","#00203FFF")) +
    theme_minimal() +
    theme(text = element_text(size = 10)) +
    ggtitle(label=variant)
  print(plot_death)
  print(kable(or_death[2:6,], caption = paste("Odds Ratio of dying for ",variant,"."),align="lllll"))
}

## Non-voc virus
cat("##Non-VOC viruses\n*Hospitalization:*\n")
freq_tab_hosp <- df_hosp %>%
  filter(Nextstrain_clade == "20A" | Nextstrain_clade=="20B"|Nextstrain_clade=="20C"|Nextstrain_clade=="20D"|Nextstrain_clade=="20G") %>%
  group_by(ORF8_ko) %>%
  summarise(n_samples = n())
print(kable(freq_tab_hosp, caption = "Non-VOC samples with hospitalization data."))
if (freq_tab_hosp[freq_tab_hosp$ORF8_ko=="yes",2]>30){
  reg_hosp <- glm(hosp ~ ORF8_ko + age_group + sex_at_birth + vaccinated, family = "binomial", data = df_hosp[df_hosp$Nextstrain_clade == "20A" | df_hosp$Nextstrain_clade=="20B"|df_hosp$Nextstrain_clade=="20C"|df_hosp$Nextstrain_clade=="20D"|df_hosp$Nextstrain_clade=="20G",])
  or_hosp <- as.tibble(questionr::odds.ratio(reg_hosp, level=(1-(0.05))), rownames = "Predictor")
  plot_hosp <- df_hosp %>%
    filter(Nextstrain_clade == "20A" | Nextstrain_clade=="20B"|Nextstrain_clade=="20C"|Nextstrain_clade=="20D"|Nextstrain_clade=="20G") %>%
    ggplot(aes(ORF8_ko, color=hosp,fill=hosp)) +
    geom_bar(position="fill",width=0.7) +
    ylab("Proportion") +
    xlab("ORF8_ko") + 
    scale_fill_manual(values =c("#ADEFD1FF","#00203FFF")) +
    scale_color_manual(values =c("#ADEFD1FF","#00203FFF")) +
    theme_minimal() +
    theme(text = element_text(size = 10)) +
    ggtitle(label="Non-VOC viruses")
  print(plot_hosp)
  print(kable(or_hosp[2:6,], caption = "Odds Ratio of Hospitalization for Non-VOC viruses.",align="lllll"))
}
cat(paste("*Death:*\n"))
freq_tab_death <- df_death %>%
  filter(Nextstrain_clade == "20A" | Nextstrain_clade=="20B"|Nextstrain_clade=="20C"|Nextstrain_clade=="20D"|Nextstrain_clade=="20G") %>%
  group_by(ORF8_ko) %>%
  summarise(n_samples = n())
print(kable(freq_tab_death, caption = "Non-VOC samples with death data."))
reg_death <- glm(died ~ ORF8_ko + age_group + sex_at_birth + vaccinated, family = "binomial", data = df_death[df_death$Nextstrain_clade == "20A" | df_death$Nextstrain_clade=="20B"|df_death$Nextstrain_clade=="20C"|df_death$Nextstrain_clade=="20D"|df_death$Nextstrain_clade=="20G",])
or_death <- as.tibble(questionr::odds.ratio(reg_death, level=(1-(0.05))), rownames = "Predictor")
plot_death <- df_death %>%
  filter(Nextstrain_clade == "20A" | Nextstrain_clade=="20B"|Nextstrain_clade=="20C"|Nextstrain_clade=="20D"|Nextstrain_clade=="20G") %>%
  ggplot(aes(ORF8_ko, color=died,fill=died)) +
  geom_bar(position="fill",width=0.7) +
  ylab("Proportion") +
  xlab("ORF8_ko") + 
  scale_fill_manual(values =c("#ADEFD1FF","#00203FFF")) +
  scale_color_manual(values =c("#ADEFD1FF","#00203FFF")) +
  theme_minimal() +
  theme(text = element_text(size = 10)) +
  ggtitle(label="Non-VOC viruses")
print(plot_death)
print(kable(or_death[2:6,], caption = "Odds Ratio of dying for Non-VOC viruses.",align="lllll"))
summary(reg_death)

reg_death_novax <- glm(died ~ ORF8_ko + age_group + sex_at_birth, family = "binomial", data = df_death[df_death$Nextstrain_clade == "20A" | df_death$Nextstrain_clade=="20B"|df_death$Nextstrain_clade=="20C"|df_death$Nextstrain_clade=="20D"|df_death$Nextstrain_clade=="20G",])
or_death_novax <- as.tibble(questionr::odds.ratio(reg_death_novax, level=(1-(0.05))), rownames = "Predictor")
kable(or_death_novax[2:5,], caption = "Odds Ratio of dying for Non-VOC viruses without vaccination status.",align="lllll")
summary(reg_death_novax)

