library(tidyverse)
library(lubridate)
library(questionr)
library(pwr)

setwd("~/Work/projects/covid/long-deletions")
## Read in SARS2 ORF8 dataset
ko <- read_tsv("results/gisaid.washington_ko_meta.tsv")
clustersAlpha <- read_tsv("nextstrain_helper/results/Alpha/clusters/clusters_ORF8.tsv")
clustersDelta <- read_tsv("nextstrain_helper/results/Delta/clusters/clusters_ORF8.tsv")
clustersOther <- read_tsv("nextstrain_helper/results/WA_other/clusters/clusters_ORF8.tsv")
wdrs <- read_csv("data/wdrs_metadata_9.2.22.csv")

## Generate cluster size
Alpha_size <- clustersAlpha %>% group_by(cluster) %>%
  add_count() %>%
  ungroup %>%
  rename(clusterSize = n) %>%
  select(strain,clusterSize)

Delta_size <- clustersDelta %>% group_by(cluster) %>%
  add_count() %>%
  ungroup %>%
  rename(clusterSize = n) %>%
  select(strain,clusterSize)

Other_size <- clustersOther %>% group_by(cluster) %>%
  add_count() %>%
  ungroup %>%
  rename(clusterSize = n) %>%
  select(strain,clusterSize)

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
  left_join(size, by = c("strain"))
            

## Plot distribution of potential predictors by ORF8_ko
# Sex at birth
df %>%
  ggplot(aes(ORF8_ko, color=sex_at_birth,fill=sex_at_birth)) +
  geom_bar(position="fill")
# Pretty even sex split in yes & no
# without omicron
df %>%
  filter(!grepl('Omicron',Nextstrain_clade)) %>%
  ggplot(aes(ORF8_ko, color=sex_at_birth,fill=sex_at_birth)) +
  geom_bar(position="fill")

# Age
df %>%
  filter(!is.na(age_group)) %>%
  ggplot(aes(ORF8_ko, color=age_group,fill=age_group)) +
  geom_bar(position="fill")
# There are younger people on average with the ko, than without. 
# May be due to variant/timing bias... But interesting...

# Age without omicron
df %>%
  filter(!is.na(age_group)) %>%
  filter(!grepl('Omicron',Nextstrain_clade)) %>%
  ggplot(aes(ORF8_ko, color=age_group,fill=age_group)) +
  geom_bar(position="fill")
## Even more age skewed without Omicron 

## without alpaha
df %>%
  filter(!is.na(age_group)) %>%
  filter(!grepl('Alpha',Nextstrain_clade)) %>%
  ggplot(aes(ORF8_ko, color=age_group,fill=age_group)) +
  geom_bar(position="fill")
## age skew mostly comes from Alpha

# Variant
df %>%
  ggplot(aes(ORF8_ko, color=Nextstrain_clade,fill=Nextstrain_clade)) +
  geom_bar(position="fill")


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
  #filter(Nextstrain_clade != "20I (Alpha, V1)") %>%
  filter(!is.na(ORF8_ko)) %>%
  ggplot(aes(days_since_vaccination,color=ORF8_ko,fill=ORF8_ko,)) +
  geom_histogram(aes(y=0.5*..density..),alpha=0.5,position='identity')
# For ORF8KO, people seem to have slightly fewer days since vacccination. 
# Should consider adding in time since shot as a continuous variable...

## Hospital
df$hosp = as_factor(df$hosp)
df$hosp = fct_relevel(df$hosp, "No", "Yes","Unknown")
df %>%
  filter(hosp!='Unknown') %>%
  filter(!grepl("Omicron",Nextstrain_clade)) %>%
  filter(!is.na(ORF8_ko)) %>%
  filter(!is.na(hosp)) %>%
  ggplot(aes(ORF8_ko, color=hosp,fill=hosp)) +
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
  filter(!grepl("Omicron",Nextstrain_clade)) %>%
  filter(!is.na(ORF8_ko)) %>%
  ggplot(aes(ORF8_ko, color=died,fill=died)) +
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
HospRegression = function(df, x) {
  df$age_group = fct_relevel(df$age_group, "0-4", "5-17", "18-44", "45-64", "65-79", "80+", "Unknown")
  df_hosp <- df %>%
    filter(hosp != "Unknown") %>%
    filter(age_group != "Unknown") %>%
    filter(!is.na(age_group)) %>%
    filter(!is.na(ORF8_ko)) %>%
    filter(sex_at_birth != 'Other') %>%
    filter(!grepl("Omicron",Nextstrain_clade)) %>%
    filter(coverage>= 0.95) %>%
    mutate(KO = ifelse(ORF8_ko == 'Yes' & clusterSize>x, 'Yes','No'))
    #mutate(KO = ifelse(ORF8_ko == 'Yes', 'Yes','No'))
  df_hosp$age_group = as.numeric(df_hosp$age_group)
  df_hosp$hosp = as.factor(df_hosp$hosp)
  
  
  reg1 <- glm(hosp ~ KO+ age_group + sex_at_birth + vaccinated, family = "binomial", data = df_hosp)
  print(summary(reg1))
  
  or_hosp = questionr::odds.ratio(reg1, level=(1-(0.05)))
  return(or_hosp)
}

for (clustSize in seq(1,20,1)) {
  print(clustSize)
  print(HospRegression(df,clustSize))
}

or_hosp_7 = HospRegression(df,7)
or_hospF = questionr::odds.ratio(hosp1,level=(1-(0.05)))

## Is time since vaccination an issue?
df$age_group = fct_relevel(df$age_group, "0-4", "5-17", "18-44", "45-64", "65-79", "80+", "Unknown")
df_vaxx <- df %>%
  filter(hosp != "Unknown") %>%
  filter(age_group != "Unknown") %>%
  filter(!is.na(age_group)) %>%
  filter(!is.na(ORF8_ko)) %>%
  filter(vaccinated=='yes') %>%
  filter(sex_at_birth != 'Other') %>%
  filter(!grepl("Omicron",Nextstrain_clade)) %>%
  filter(coverage>= 0.95)
df_vaxx$age_group = as.numeric(df_vaxx$age_group)
df_vaxx$hosp = as.factor(df_vaxx$hosp)
  
vax1 <- glm(hosp ~ ORF8_ko + age_group + sex_at_birth + days_since_vaccination, family = "binomial", data = df_vaxx)
print(summary(vax1))
vax2 <- glm(hosp ~ ORF8_ko + age_group + sex_at_birth, family = "binomial", data = df_vaxx)
summary(vax2)
vax3 <- glm(hosp ~ ORF8_ko*days_since_vaccination + age_group + sex_at_birth, family = "binomial", data = df_vaxx)
summary(vax3)
vax4 <- glm(hosp ~ days_since_vaccination + age_group + sex_at_birth, family = "binomial", data = df_vaxx)
summary(vax4)
vax5 <- glm(hosp ~ days_since_vaccination*age_group + sex_at_birth, family = "binomial", data = df_vaxx)
summary(vax5)
vax6 <- glm(hosp ~ age_group + sex_at_birth*days_since_vaccination, family = "binomial", data = df_vaxx)
summary(vax6)
## ORF8KO doesn't seem to have an effect on clinical outcomes in vaccinated people.
## instead, i'm observing this weird effect in which the more days out you are
## from vaccination, the less likely you are to get sick...

df_hosp <- df %>%
  filter(hosp != "Unknown") %>%
  filter(age_group != "Unknown") %>%
  filter(!is.na(age_group)) %>%
  filter(!is.na(ORF8_ko)) %>%
  filter(sex_at_birth != 'Other') %>%
  filter(!grepl("Omicron",Nextstrain_clade)) %>%
  filter(coverage>= 0.95)
df_hosp$age_group = as.numeric(df_hosp$age_group)
df_hosp$hosp = as.factor(df_hosp$hosp)

hosp1 <- glm(hosp ~ ORF8_ko + age_group + sex_at_birth + vaccinated, family = "binomial", data = df_hosp)
summary(hosp1)
hosp2 <- glm(hosp ~ ORF8_ko*vaccinated + age_group + sex_at_birth , family = "binomial", data = df_hosp)
summary(hosp2)
hosp3 <- glm(hosp ~ ORF8_ko + age_group*vaccinated + sex_at_birth , family = "binomial", data = df_hosp)
summary(hosp3)
hosp4 <- glm(hosp ~ ORF8_ko + age_group +vaccinated*sex_at_birth , family = "binomial", data = df_hosp)
summary(hosp4)
hosp5 <- glm(hosp ~ ORF8_ko*age_group + sex_at_birth + vaccinated, family = "binomial", data = df_hosp)
summary(hosp5)
hosp6 <- glm(hosp ~ ORF8_ko*sex_at_birth + age_group +  vaccinated, family = "binomial", data = df_hosp)
summary(hosp6)
hosp7 <- glm(hosp ~ ORF8_ko*age_group*vaccinated + sex_at_birth, family = "binomial", data = df_hosp)
summary(hosp7)
hosp8 <- glm(hosp ~ ORF8_ko:age_group + ORF8_ko + age_group + ORF8_ko:vaccinated + vaccinated + sex_at_birth, family = "binomial", data = df_hosp)
summary(hosp8)

## Effect of ORF8KO & vaccines seems to be stronger in younger people...
## Effect of ORF8KO only in unvaccinated people
## Timing of vaccine rollout?
## Variant timing??
df_noVax <- df %>%
  filter(hosp != "Unknown") %>%
  filter(age_group != "Unknown") %>%
  filter(!is.na(age_group)) %>%
  filter(!is.na(ORF8_ko)) %>%
  filter(sex_at_birth != 'Other') %>%
  filter(!grepl("Omicron",Nextstrain_clade)) %>%
  filter(coverage>= 0.95) %>%
  filter(vaccinated=='no')
df_noVax$age_group = as.numeric(df_noVax$age_group)
df_noVax$hosp = as.factor(df_noVax$hosp)

noVax1 <- glm(hosp ~ ORF8_ko +age_group + sex_at_birth, family = "binomial", data = df_noVax)
summary(noVax1)
noVax2 <- glm(hosp ~ ORF8_ko*age_group + sex_at_birth, family = "binomial", data = df_noVax)
summary(noVax2)
## Yeah biggest impact in unvaccinated


# Death
DeathRegression = function(df,x) {
  df_death <- df %>%
  filter(age_group != "Unknown") %>%
  filter(!is.na(age_group)) %>%
  filter(!is.na(ORF8_ko)) %>%
  filter(!is.na(died)) %>% 
  filter(sex_at_birth != 'Other') %>%
  filter(!grepl("Omicron",Nextstrain_clade)) %>%
  filter(coverage>= 0.95) %>%
  filter(!is.na(sex_at_birth)) %>%
  #mutate(KO = ifelse(ORF8_ko == 'Yes' & clusterSize>x, 'Yes','No'))
  mutate(KO = ifelse(ORF8_ko == 'Yes' & clusterSize>x, 'Yes','No'))
  df_death$age_group = as.numeric(df_death$age_group)
  df_death$died = fct_relevel(df_death$died, "No", "Yes")
  
  reg2 <- glm(died ~ KO + age_group + sex_at_birth + vaccinated, family = "binomial", data = df_death)
  print(summary(reg2))
  
  or_death = questionr::odds.ratio(reg2, level=(1-(0.05)))
  return(or_death)
}  


for (clustSize in seq(1,20,1)) {
  print(clustSize)
  print(DeathRegression(df,clustSize))
}

df_death <- df %>%
  filter(age_group != "Unknown") %>%
  filter(!is.na(age_group)) %>%
  filter(!is.na(ORF8_ko)) %>%
  filter(!is.na(died)) %>% 
  filter(sex_at_birth != 'Other') %>%
  filter(!grepl("Omicron",Nextstrain_clade)) %>%
  filter(coverage>= 0.95) %>%
  filter(!is.na(sex_at_birth)) %>%
  filter(!is.na(vaccinated))
df_death$age_group = as.numeric(df_death$age_group)
df_death$died = fct_relevel(df_death$died, "No", "Yes")


death1 <- glm(died ~ ORF8_ko + age_group + sex_at_birth + vaccinated, family = "binomial", data = df_death)
summary(death1)
death2 <- glm(died ~ ORF8_ko + age_group+ vaccinated, family = "binomial", data = df_death)
summary(death2)

summary(death1)
or_death <- DeathRegression(df,1)
jtools::summ(hosp7)

### Do we even have power to detect an effect of death?
df_death %>%
  group_by(ORF8_ko, died) %>%
  count()

n1 = 8203 + 107
n2 = 40362+796
p2 = 796/n2
p1 = p2*0.86
pwr.2p2n.test(h=ES.h(p1,p2),n1 =n1, n2 =n2, alternative='less',sig.level=0.05)
#pwr.f2.test(u=4,v=49457,f2=0.13)

or_death <- cbind(Variable = rownames(or_death), or_death)

  
or_hospF <- cbind(Variable = rownames(or_hospF), or_hospF)
or_hospF <- cbind(ko = c('no','yes','no','no','no'), or_hospF)

or_death <- cbind(ko = c('no','yes','no', 'no','no'), or_death)
colnames(or_hospF)[4] ="minCI"
colnames(or_hospF)[5] ="maxCI"
colnames(or_death)[4] ="minCI"
colnames(or_death)[5] ="maxCI"


or_hospF %>%
  filter(Variable!="(Intercept)") %>%
  ggplot(aes(y=Variable,x=OR)) +
  geom_vline(xintercept=1,linetype='dashed') +
  geom_linerange(aes(xmin=minCI,xmax=maxCI)) +
  geom_point(size=5,aes(color=ko,fill=ko))+
  scale_x_log10() +
  theme_minimal() +
  scale_y_discrete(labels=c('Age group', 'ORF8 KO','Sex: Male', 'Vaccinated')) +
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
  scale_y_discrete(labels=c('Age group', 'ORF8 KO','Sex: Male','Vaccinated')) +
  ylab('') +
  xlab('Odds ratio for death') +
  scale_color_manual(values=c('grey','#ce4257')) +
  scale_fill_manual(values=c('grey','#ce4257')) +
  theme(legend.position = 'none')
ggsave('figs/or_death.pdf',dpi=300,height=3,width=4)
                       








