library(tidyverse)
library(lubridate)
library(questionr)
library(pwr)
library(ggsignif)
library(ggpubr)
library(patchwork)
library(table1)

setwd("~/Work/projects/covid/long-deletions")
## Read in SARS2 ORF8 dataset
df <- read_tsv("wa_results/clinical_combined.tsv")

## for demo use:
#df <- read_tsv('data/clinical_example.tsv')



df$age_group = fct_relevel(df$age_group, "0-4", "5-17", "18-44", "45-64", "65-79", "80+")



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
df %>%
  #filter(Nextstrain_clade != "20I (Alpha, V1)") %>%
  filter(!is.na(ORF8_ko)) %>%
  ggplot(aes(days_since_vaccination,color=ORF8_ko,fill=ORF8_ko,)) +
  geom_histogram(aes(y=0.5*..density..),alpha=0.5,position='identity')
# For ORF8KO, people seem to have slightly fewer days since vacccination.
# Should consider adding in time since shot as a continuous variable...

## Hospital
df$hosp = factor(df$hosp,levels=c("No","Yes"))
#df$hosp = fct_relevel(df$hosp, "No", "Yes")
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

hosp.2way <- df %>%
  filter(hosp!='Unknown') %>%
  filter(!grepl("Omicron",Nextstrain_clade)) %>%
  filter(!is.na(ORF8_ko)) %>%
  filter(!is.na(hosp)) %>%
  with(table(hosp, ORF8_ko)) %>% prop.table(margin=2) %>%
  as_tibble()

df %>%
  filter(hosp!='Unknown') %>%
  filter(!grepl("Omicron",Nextstrain_clade)) %>%
  filter(!is.na(ORF8_ko)) %>%
  filter(!is.na(hosp)) %>%
  with(table(hosp, ORF8_ko))

hosp.2way.mat <- df %>%
  filter(!grepl("Omicron",Nextstrain_clade)) %>%
  filter(!is.na(ORF8_ko)) %>%
  filter(!is.na(hosp)) %>%
  filter(hosp!='Unknown') %>%
  with(table(hosp, ORF8_ko))


## Death

df$died = as_factor(df$died)
df$died = fct_relevel(df$died, "No", "Yes")

df %>%
  filter(!is.na(ORF8_ko)) %>%
  filter(!is.na(died)) %>%
  filter(!grepl("Omicron",Nextstrain_clade)) %>%
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


died.2way <- df %>%
  filter(!is.na(died)) %>%
  filter(!grepl("Omicron",Nextstrain_clade)) %>%
  filter(!is.na(ORF8_ko)) %>%
  with(table(died, ORF8_ko)) %>% prop.table(margin=2) %>%
  as_tibble()

df %>%
  filter(!is.na(died)) %>%
  filter(!grepl("Omicron",Nextstrain_clade)) %>%
  filter(!is.na(ORF8_ko)) %>%
  with(table(died, ORF8_ko))

died.2way.mat <- df %>%
  filter(!is.na(died)) %>%
  filter(!grepl("Omicron",Nextstrain_clade)) %>%
  filter(!is.na(ORF8_ko)) %>%
  with(table(died, ORF8_ko))

deathp <- chisq.test(died.2way.mat)$p.value
hospp <- chisq.test(hosp.2way.mat)$p.value
fisher.test(died.2way.mat)$p.value
fisher.test(hosp.2way.mat)$p.value

clin.prop <- hosp.2way %>%
  full_join(died.2way) %>%
  pivot_longer(cols=c('hosp','died'),names_to='outcome') %>%
  mutate(outcome = as_factor(outcome)) %>%
  mutate(outcome = fct_relevel(outcome,c('hosp','died'))) %>%
  filter(!is.na(value)) %>%
  filter(value!='No') %>%
  ggplot(aes(x=ORF8_ko,y=n,fill=outcome,color=outcome)) +
  geom_bar(stat='identity',position=position_dodge()) +
  ylab('Proportion') +
  xlab('ORF8 knockout') +
  scale_fill_manual(labels = c('Hospitalized','Died'),values = c('#457b9d','#14213d'),name='Outcome:') +
  scale_color_manual(labels = c('Hospitalized','Died'),values = c('#457b9d','#14213d'),name='Outcome:') +
  theme_minimal() +
  geom_signif(
    y_position = c(0.069, 0.087), xmin = c(1.25, 0.75), xmax = c(2.25, 1.75),
    annotation = c(signif(deathp,digits=2), signif(hospp,digits=0)),color='black'
  ) +
  theme(legend.position='bottom') +
  ylim(0,0.09)
clin.prop
ggsave('figs/fig5/proportions_clinical.pdf',dpi=300,height=3.5,width=3)

hosp.2way %>%
  full_join(died.2way) %>%
  pivot_longer(cols=c('hosp','died'),names_to='outcome') %>%
  mutate(outcome = as_factor(outcome)) %>%
  mutate(outcome = fct_relevel(outcome,c('hosp','died'))) %>%
  filter(!is.na(value)) %>%
  filter(value!='No') %>%
  rename(Proportion = n) %>%
  select(ORF8_ko,Proportion,outcome) %>%
  write_tsv('figs/fig5/5A_SourceData.tsv')


df_agg <- df %>%
  filter(age_group != "Unknown") %>%
  filter(!is.na(age_group)) %>%
  filter(!is.na(ORF8_ko)) %>%
  filter(sex_at_birth != 'Other') %>%
  filter(!grepl("Omicron",Nextstrain_clade)) %>%
  filter(coverage>= 0.95) %>%
  filter(!is.na(sex_at_birth)) %>%
  mutate(hosp = if_else(is.na(hosp),'Unknown',hosp)) %>%
  filter(!is.na(vaccinated)) %>%
  mutate(VOC = if_else(grepl("Alpha",Nextstrain_clade)|grepl("Delta",Nextstrain_clade)|grepl("Beta",Nextstrain_clade)|grepl("Gamma",Nextstrain_clade),"Yes","No"))

df_agg$died = fct_relevel(df_agg$died, "No", "Yes")

df_agg %>%
  select(ORF8_ko,sex_at_birth,vaccinated, VOC,hosp,died) %>%
  write_tsv('data/clinical_subset.tsv')

## Table
no_omicron <- df_agg %>%
  filter(!grepl("Omicron",Nextstrain_clade)) %>%
  mutate(vaccinated = ifelse(vaccinated=='no','No','Yes')) %>%
  mutate(ORF8_ko = ifelse(ORF8_ko == 'Yes', 'ORF8 knockout','ORF8 intact')) %>%
  dplyr::select(died,hosp,ORF8_ko,vaccinated,age_group,sex_at_birth,VOC,wkyr)

label(no_omicron$ORF8_ko) = 'ORF8 knockout?'
label(no_omicron$died) = 'Died?'
label(no_omicron$hosp) = 'Hospitalized?'
label(no_omicron$vaccinated) = 'Vaccinated'
label(no_omicron$sex_at_birth) = 'Sex assigned at birth'
label(no_omicron$VOC) = 'Variant of concern?'
no_omicron$age_group = factor(no_omicron$age_group, levels=c("0-4", "5-17", "18-44", "45-64", "65-79", "80+"))
label(no_omicron$age_group) = 'Age group'
label(no_omicron$wkyr) = 'Year + Epiweek'
table1(~hosp + died + VOC + age_group + sex_at_birth + vaccinated | ORF8_ko, data=no_omicron,topclass="Rtable1-grid")



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
    mutate(KO = ifelse(ORF8_ko == 'Yes' & clusterSize>x, 'Yes','No')) %>%
    mutate(VOC = if_else(grepl("Alpha",Nextstrain_clade)|grepl("Delta",Nextstrain_clade)|grepl("Beta",Nextstrain_clade)|grepl("Gamma",Nextstrain_clade),"Yes","No"))
    #mutate(KO = ifelse(ORF8_ko == 'Yes', 'Yes','No'))
  df_hosp$age_group = as.numeric(df_hosp$age_group)
  df_hosp$hosp = as.factor(df_hosp$hosp)


  reg1 <- glm(hosp ~ KO+ age_group + sex_at_birth + vaccinated + VOC + wkyr, family = "binomial", data = df_hosp)
  #print(summary(reg1))

  or_hosp = cbind(OR = exp(reg1$coefficients[2:6]), exp(confint(reg1,parm=c("KOYes","age_group","sex_at_birthMale","vaccinatedyes","VOCYes"))))

  results = list(or_hosp,reg1$aic)
  return(results)
}

## Is time since vaccination an issue?
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
  filter(coverage>= 0.95) %>%
  mutate(VOC = if_else(grepl("Alpha",Nextstrain_clade)|grepl("Delta",Nextstrain_clade)|grepl("Beta",Nextstrain_clade)|grepl("Gamma",Nextstrain_clade),"Yes","No"))
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
hosp9 <- glm(hosp ~ ORF8_ko + age_group + sex_at_birth + vaccinated + VOC, family = binomial, data = df_hosp)
summary(hosp9)
hosp10 <- glm(hosp ~ ORF8_ko*VOC + age_group + sex_at_birth + vaccinated, family = binomial, data = df_hosp)
summary(hosp10)

hosp11 <- glm(hosp ~ ORF8_ko + age_group + sex_at_birth + vaccinated + VOC + wkyr,family = binomial, data = df_hosp)
summary(hosp11)

hosp12 <- glm(hosp ~ ORF8_ko + age_group + sex_at_birth + vaccinated + VOC + moyr,family = binomial, data = df_hosp)
summary(hosp12)

hosp13 <- glm(hosp ~ ORF8_ko + age_group + sex_at_birth + vaccinated + wkyr + VOC,family = binomial, data = df_hosp[df_hosp$Nextstrain_clade!='20I (Alpha, V1)',])
summary(hosp13)
## A lot of signal is coming from Alpha

hosp14 <- glm(hosp ~ age_group + sex_at_birth + vaccinated + VOC + wkyr,family = binomial, data = df_hosp)
summary(hosp14)

hosp15 <- glm(hosp ~ ORF8_ko + sex_at_birth + vaccinated + VOC, family=binomial, data = df_hosp)
summary(hosp15)

df_hosp_spec <- df_hosp %>% mutate(age_group = ifelse(age_group==1|age_group==2|age_group==3, 3, age_group))
hosp16 <- glm(hosp ~ ORF8_ko + age_group + sex_at_birth + vaccinated + VOC + wkyr,family = binomial, data = df_hosp_spec)
summary(hosp16)

hosp17 <- glm(hosp ~ ORF8_ko + age_group + VOC,family = binomial, data = df_hosp_spec)
summary(hosp17)

hosp18 <- glm(hosp ~ ORF8_ko + age_group + VOC,family = binomial, data = df_hosp)
summary(hosp18)

exp(hosp17$coefficients)
exp(hosp15$coefficients)
exp(hosp11$coefficients)
exp(hosp18$coefficients)

or_hospF = cbind(OR = exp(hosp11$coefficients[2:6]), exp(confint(hosp11,parm=c("ORF8_koYes","age_group","sex_at_birthMale","vaccinatedyes","VOCYes"))))

## power calculations for hospital
r2A = pseudoR2(hosp14)
r2AB = pseudoR2(hosp11)
f2 = (r2AB - r2A)/(1-r2AB)
v=25460
u = 25531-v-1
pwr.f2.test(u=u,v=v,f2=f2)


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
  mutate(KO = ifelse(ORF8_ko == 'Yes' & clusterSize>x, 'Yes','No')) %>%
    mutate(VOC = if_else(grepl("Alpha",Nextstrain_clade)|grepl("Delta",Nextstrain_clade)|grepl("Beta",Nextstrain_clade)|grepl("Gamma",Nextstrain_clade),"Yes","No"))
  df_death$age_group = as.numeric(df_death$age_group)
  df_death$died = fct_relevel(df_death$died, "No", "Yes")

  reg2 <- glm(died ~ KO + age_group + sex_at_birth + vaccinated + VOC + wkyr, family = "binomial", data = df_death)
  #print(summary(reg2))

  or_death = cbind(OR = exp(reg2$coefficients[2:6]), exp(confint(reg2,parm=c("KOYes","age_group","sex_at_birthMale","vaccinatedyes","VOCYes"))))

  results = list(or_death, reg2$aic)
  return(results)
}

HospRegression(df,1)

for (clustSize in seq(1,50,1)) {
  if (exists("effect_size")) {
    reg1 <- DeathRegression(df,clustSize)
    or1 <- cbind(variable = rownames(reg1[[1]]), reg1[[1]],clusterSize = clustSize,aic = reg1[[2]],regression='death')
    effect_size <- rbind(effect_size,or1)
  } else {
    reg1 <- DeathRegression(df,clustSize)
    effect_size <- cbind(variable = rownames(reg1[[1]]), reg1[[1]], clusterSize = clustSize,aic = reg1[[2]], regression ='death')
  }
  reg2 <- HospRegression(df,clustSize)
  or2 <- cbind(variable = rownames(reg2[[1]]), reg2[[1]], clusterSize = clustSize,aic = reg2[[2]],regression='hospitalization')
  effect_size <- rbind(effect_size,or2)
}

effect <- as_tibble(effect_size)
effect$OR = as.numeric(effect$OR)
effect$`2.5 %` = as.numeric(effect$`2.5 %`)
effect$`97.5 %` = as.numeric(effect$`97.5 %`)
effect$clusterSize = as.numeric(effect$clusterSize)
effect$aic = as.numeric(effect$aic)

effect %>% 
  filter(variable=='KOYes') %>%
  select(-variable) %>%
  write_tsv('figs/supplemental/S9_SourceData.tsv')

odds_cluster <- effect %>%
  filter(variable=='KOYes') %>%
  ggplot(aes(x=clusterSize,y=OR,fill=regression,color=regression,)) +
  geom_line() +
  geom_ribbon(aes(ymin=`2.5 %`, ymax=`97.5 %`),alpha=0.3,color=NA) +
  scale_fill_manual(labels = c('Death','Hospitalization'),values = c('#14213d','#457b9d')) +
  scale_color_manual(labels = c('Death','Hospitalization'),values = c('#14213d','#457b9d')) +
  ylab('Odds ratio for ORF8 knockout') +
  xlab('') +
  theme_minimal() +
  theme(axis.ticks.x = element_blank(), axis.text.x=element_blank(),
        legend.position = 'top',legend.title=element_blank()) +
  xlim(0,50)

odds_cluster

scale <- 0.585
aic_cluster <-effect %>%
  filter(variable=='KOYes') %>%
  dplyr::select(clusterSize,regression,aic) %>%
  pivot_wider(names_from = regression, values_from = aic) %>%
  ggplot(aes(x=clusterSize)) +
  geom_line(aes(y=hospitalization),color='#457b9d') +
  geom_line(aes(y=death/scale),color='#14213d') +
  #scale_fill_manual(labels = c('Death','Hospitalization'),values = c('#14213d','#457b9d')) +
  #scale_color_manual(labels = c('Death','Hospitalization'),values = c('#14213d','#457b9d')) +
  xlab('Cluster size required \nto be defined as knockout') +
  theme_minimal() +
  scale_y_continuous(sec.axis = sec_axis(~.*scale, name="AIC for death model")) +
  ylab('AIC for hospitalization model') +
  theme(
    axis.title.y = element_text(color = "#457b9d"),
    axis.title.y.right = element_text(color = "#14213d")
  ) +
  xlim(0,50)


aic_cluster


odds_cluster + aic_cluster + plot_layout(nrow=2) + plot_annotation(tag_levels = 'A') & scale_x_continuous(limits = c(0, 50))

ggsave('figs/supplemental/OddsRatio_clustersize_time.pdf',dpi=300,height=6,width=4)
ggsave('figs/supplemental/OddsRatio_clustersize_time.jpg',dpi=300,height=6,width=4,bg='white')

df_death <- df %>%
  filter(age_group != "Unknown") %>%
  filter(!is.na(age_group)) %>%
  filter(!is.na(ORF8_ko)) %>%
  filter(!is.na(died)) %>%
  filter(sex_at_birth != 'Other') %>%
  filter(!grepl("Omicron",Nextstrain_clade)) %>%
  filter(coverage>= 0.95) %>%
  filter(!is.na(sex_at_birth)) %>%
  filter(!is.na(vaccinated)) %>%
  mutate(VOC = if_else(grepl("Alpha",Nextstrain_clade)|grepl("Delta",Nextstrain_clade)|grepl("Beta",Nextstrain_clade)|grepl("Gamma",Nextstrain_clade),"Yes","NO"))

df_death$age_group = as.numeric(df_death$age_group)
df_death$died = fct_relevel(df_death$died, "No", "Yes")








death1 <- glm(died ~ ORF8_ko + age_group + sex_at_birth + vaccinated, family = "binomial", data = df_death)
summary(death1)
death2 <- glm(died ~ ORF8_ko + age_group+ vaccinated, family = "binomial", data = df_death)
summary(death2)
death3 <- glm(died ~ ORF8_ko + age_group + sex_at_birth + vaccinated + VOC, family = "binomial", data = df_death)
summary(death3)

death4 <- glm(died ~ ORF8_ko + age_group + sex_at_birth + vaccinated + VOC + wkyr, family = "binomial", data = df_death)
summary(death4)

death5<- glm(died ~ ORF8_ko + age_group + sex_at_birth + vaccinated + wkyr,family = binomial, data = df_death[df_death$Nextstrain_clade!='20I (Alpha, V1)',])
summary(death5)

death6 <- glm(died ~ age_group + sex_at_birth + vaccinated + VOC + wkyr, family = "binomial", data = df_death)
summary(death4)

death7 <- glm(died ~ ORF8_ko + sex_at_birth + vaccinated + VOC, family=binomial, data = df_death)
summary(death7)

df_death_spec <- df_death %>% mutate(age_group = ifelse(age_group==1|age_group==2|age_group==3, 3, age_group))
death8 <- glm(died ~ ORF8_ko + age_group + sex_at_birth + vaccinated + VOC + wkyr, family = "binomial", data = df_death_spec)
summary(death8)

death9 <- glm(died ~ ORF8_ko + age_group + VOC, family = "binomial", data = df_death_spec)
summary(death9)

## Lose effect if you eliminate Alpha, but not enough power?

cbind(OR = exp(death4$coefficients[2:6]), exp(confint(death4,parm=c("ORF8_koYes","age_group","sex_at_birthMale","vaccinatedyes","VOCYes"))))


or_death = cbind(OR = exp(death4$coefficients[2:6]), exp(confint(death4,parm=c("ORF8_koYes","age_group","sex_at_birthMale","vaccinatedyes","VOCYes"))))

### Do we even have power to detect an effect of death?
pseudoR2 = function(model){
  return(with(summary(model), 1 - deviance/null.deviance))
}


df_death %>%
  group_by(ORF8_ko, died) %>%
  count()

r2A = pseudoR2(death6)
r2AB = pseudoR2(death4)
f2 = (r2AB - r2A)/(1-r2AB)
v=49385
u = 49468-v-1
pwr.f2.test(u=u,v=v,f2=f2)

n1 = 8203 + 107
n2 = 40362+796
p2 = 796/n2
p1 = p2*0.73
pwr.2p2n.test(h=ES.h(p1,p2),n1 =n1, n2 =n2, alternative='less',sig.level=0.05)
#pwr.f2.test(u=4,v=49457,f2=0.13)

## Potential omicron power
### Do we even have power to detect an effect of death?

n1 = 4820+323
n2 = 19627 + 1629
p2 = 1629/n2
p1 = p2*0.82
pwr.2p2n.test(h=ES.h(p1,p2),n1 =n1, n2 =n2, alternative='less',sig.level=0.05)
#pwr.f2.test(u=4,v=49457,f2=0.13)


## Plot odds ratios
or_death <- cbind(Variable = rownames(or_death), or_death)
or_hospF <- cbind(Variable = rownames(or_hospF), or_hospF)
colnames(or_hospF)[3] ="minCI"
colnames(or_hospF)[4] ="maxCI"
colnames(or_death)[3] ="minCI"
colnames(or_death)[4] ="maxCI"
or_death <- cbind(or_death, regression='death')
or_hospF <- cbind(or_hospF, regression='hospitalization')
all_odds <- rbind(or_death, or_hospF) %>%
  as_tibble()

all_odds$OR = as.numeric(all_odds$OR)
all_odds$minCI = as.numeric(all_odds$minCI)
all_odds$maxCI = as.numeric(all_odds$maxCI)

all_odds %>%
  write_tsv('figs/fig5/5B_SourceData.tsv')

odds_both <- all_odds %>%
  ggplot(aes(y=Variable,x=OR)) +
  geom_vline(xintercept=1,linetype='dashed') +
  geom_pointrange(size=0.7,aes(xmin=minCI, xmax=maxCI, color=regression,fill=regression),position=position_dodge(width=0.6))+
  scale_x_log10() +
  theme_minimal() +
  scale_y_discrete(labels=c('Age group increase', 'ORF8 knockout','Sex: Male', 'Vaccination','Variant of Concern')) +
  scale_fill_manual(labels = c('Death','Hospitalization'),values = c('#14213d','#457b9d')) +
  scale_color_manual(labels = c('Death','Hospitalization'),values = c('#14213d','#457b9d')) +
  ylab('') +
  xlab('Odds ratio') +
  theme(axis.text.y=element_text(size=11,color='black'))
odds_both
ggsave('figs/fig5/OddsRatio_time.pdf',dpi=300,height=4,width=6)
ggsave('figs/fig5/OddsRatio_time.jpg',dpi=300,height=4,width=6)

odds_both_no_legend <- odds_both + theme(legend.position='none')

ggarrange(clin.prop,odds_both_no_legend,labels=c('A','B'),widths=c(0.5,1))#,common.legend=TRUE,legend='bottom')
ggsave('figs/fig5/fig5_time.pdf',dpi=300,height=3,width=6.5)
ggsave('figs/fig5/fig5_time.jpg',dpi=300,height=3,width=6.5,bg='white')


### 19B specific
next19b_hosp = next19b %>%
  filter(hosp != "Unknown") %>%
  #filter(age_group != "Unknown") %>%
  #filter(!is.na(age_group)) %>%
  filter(!is.na(ORF8_ko)) %>%
  #filter(sex_at_birth != 'Other') %>%
  filter(coverage>= 0.95)
next19b_hosp$age_group = as.numeric(next19b_hosp$age_group)
next19b_hosp$hosp = as.factor(next19b_hosp$hosp)


next19b_death <- next19b %>%
  #filter(age_group != "Unknown") %>%
  #filter(!is.na(age_group)) %>%
  filter(!is.na(ORF8_ko)) %>%
  filter(!is.na(died)) %>%
  #filter(sex_at_birth != 'Other') %>%
  filter(coverage>= 0.95)# %>%
  #filter(!is.na(sex_at_birth))
next19b_death$age_group = as.numeric(next19b_death$age_group)
next19b_death$died = fct_relevel(next19b_death$died, "No", "Yes")
