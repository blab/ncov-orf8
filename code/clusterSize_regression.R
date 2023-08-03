## Runs poisson regression on cluster size

# Load packages
library(tidyverse)
library(lubridate)
library(MASS)
library(ggplot2)


#  Change working directory
setwd("~/Work/projects/covid/long-deletions")

# Read in data
load_data = function(path,ober) {
  df = read_tsv(path) %>% filter(mut_type!='undoStop')
  df$date_observed = ymd(df$date_observed)
  df$time = interval(df$date_observed, ymd(20230501)) %>% as.numeric("years")
  df <- df %>% left_join(ober, by=c("aa_mutations"="mutation")) %>%
    mutate(`Δ log R` = if_else(is.na(`Δ log R`),0,`Δ log R`)) %>%
    mutate(`Δ log R 95% ci upper` = if_else(is.na(`Δ log R 95% ci upper`),0,`Δ log R 95% ci upper`)) %>%
    mutate(`Δ log R 95% ci lower` = if_else(is.na(`Δ log R 95% ci lower`),0,`Δ log R 95% ci lower`)) %>%
    mutate(n_descendants = leaf_count - 1) %>%
    mutate(fitness = if_else(`Δ log R 95% ci lower`>0,'increased','same')) %>%
    mutate(mut_type_obermeyer = ifelse(fitness=='increased','missense_increased',mut_type)) %>%
    mutate(mut_type_obermeyer = ifelse(mut_type_obermeyer=='missense','missense_other',mut_type_obermeyer))
  df$mut_type_obermeyer = as_factor(df$mut_type_obermeyer)
  df$mut_type_obermeyer = fct_relevel(df$mut_type_obermeyer,c('synonymous','missense_increased','missense_other','nonsense'))
  df$mut_type = as_factor(df$mut_type)
  df$mut_type = fct_relevel(df$mut_type, c('synonymous','missense','nonsense'))
  return(df)
}
obermeyer = read_tsv('data/obermeyer_mutations.tsv')
orf8 = load_data('usher/trimmed/clades_nested/ORF8_clades.tsv',obermeyer)
orf1a = load_data('usher/trimmed/clades_nested/ORF1a_clades.tsv',obermeyer)
spike = load_data('usher/trimmed/clades_nested/S_clades.tsv',obermeyer)
orf1b = load_data('usher/trimmed/clades_nested/ORF1b_clades.tsv',obermeyer)
orf3a = load_data('usher/trimmed/clades_nested/ORF3a_clades.tsv',obermeyer)
e = load_data('usher/trimmed/clades_nested/E_clades.tsv',obermeyer)
m = load_data('usher/trimmed/clades_nested/M_clades.tsv',obermeyer)
orf6 = load_data('usher/trimmed/clades_nested/ORF6_clades.tsv',obermeyer)
orf7a = load_data('usher/trimmed/clades_nested/ORF7a_clades.tsv',obermeyer)
orf7b = load_data('usher/trimmed/clades_nested/ORF7b_clades.tsv',obermeyer)
n = load_data('usher/trimmed/clades_nested/N_clades.tsv',obermeyer)
orf9b = load_data('usher/trimmed/clades_nested/ORF9b_clades.tsv',obermeyer)


# Run regression
run_pois = function(df) {
  df %>%
    mutate(clusterSize = leaf_count - 1) -> filt
  reg <- glm(clusterSize ~ mut_type, family = poisson(link="log"),data=filt,offset=log(time))
  print(summary(reg))
  print(exp(cbind(Odds=coef(reg), confint(reg))))
  return(reg)
}

run_qpois = function(df) {
  df %>%
    filter(mut_type != 'undoStop') %>%
    mutate(clusterSize = leaf_count - 1) -> filt
  reg <- glm(clusterSize ~ mut_type, family = quasipoisson(link="log"),data=filt,offset=log(time))
  print(summary(reg))
  print(exp(cbind(Odds=coef(reg), confint(reg))))
  return(reg)
}

orf8_pois = run_pois(orf8)
orf8_qpois = run_qpois(orf8)

nb_reg = function(df,maxit) {
  fit = glm.nb(n_descendants ~ mut_type + offset(log(time)),data=df,maxit=maxit)
  print(summary(fit))
  return(fit)
}

nb_reg_obermeyer = function

coefficients =  function(fit,gene){
  df = exp(cbind(Fold_change=coef(fit),confint(fit)))
  df <- cbind(Variable = rownames(df), df)
  df <- cbind(df,Gene=gene)
  return(as_tibble(df))
}

orf8_nb = nb_reg(orf8,50)
orf8_coef = coefficients(orf8_nb,'ORF8')
orf1a_nb = nb_reg(orf1a,100)
orf1a_coef = coefficients(orf1a_nb,'ORF1a')
spike_nb = nb_reg(spike,200)
spike_coef = coefficients(spike_nb,'Spike')  
orf1b_nb = nb_reg(orf1b,100)
orf1b_coef = coefficients(orf1b_nb,'ORF1b')
orf3a_nb = nb_reg(orf3a,100)
orf3a_coef = coefficients(orf3a_nb,'ORF3a')  
e_nb = nb_reg(e,200)
e_coef = coefficients(e_nb,'E')  ## THIS DOESN't WORK
m_nb = nb_reg(m,300)
m_coef = coefficients(m_nb,'M')  
orf6_nb = nb_reg(orf6,300)
orf6_coef = coefficients(orf6_nb,'ORF6')
orf7a_nb = nb_reg(orf7a,100)
orf7a_coef = coefficients(orf7a_nb,'ORF7a')
orf7b_nb = nb_reg(orf7b,300)
orf7b_coef = coefficients(orf7b_nb,'ORF7b')  
n_nb = nb_reg(n,200)
n_coef = coefficients(n_nb,'N')  
orf9b_nb = nb_reg(orf9b,200)
orf9b_coef = coefficients(orf9b_nb,'ORF9b')  

results = bind_rows(orf1a_coef,orf1b_coef,spike_coef,orf3a_coef,m_coef,orf6_coef,orf7a_coef,orf7b_coef,orf8_coef,n_coef,orf9b_coef)

orf1a_obermeyer_ceof = coefficients(orf1a_obermeyer_nb,'ORF1a')

results_small = bind_rows(orf1a_obermeyer_ceof,spike_coef,orf8_coef,orf1a_coef)

results_small$Fold_change = as.numeric(results_small$Fold_change)
results_small$`2.5 %` = as.numeric(results_small$`2.5 %`)
results_small$`97.5 %` = as.numeric(results_small$`97.5 %`)
results_small$Gene = as_factor(results_small$Gene)
results_small$Gene = fct_relevel(results_small$Gene,c('ORF1a','Spike','ORF8'))
mapping = c("mut_type_obermeyermissense_increased"="Missense: Increased fitness\n(Obermeyer et al)", "mut_type_obermeyermissense_other"="Missense: Other",'mut_type_obermeyernonsense'="Nonsense","mut_typemissense"="Missense","mut_typenonsense"="Nonsense")



## Palette
pal1 = c('#e71d36','#2EC4B6','#ff9f1c','#9E641B','#FFEBD2','#A2E1D8','#266962')
colorblindcheck::palette_check(pal1, plot = TRUE)


results_small %>%
  filter(Variable!='(Intercept)') %>%
  filter(Gene == 'Spike'|Gene =='ORF1a'|Gene =='ORF8') %>%
  mutate(minci = replace_na(`2.5 %`,0.001)) %>%
  mutate(maxci = replace_na(`97.5 %`,9.99)) %>%
  mutate(Mutation_type = mapping[Variable])%>% 
  mutate(Mutation_type = fct_relevel(Mutation_type, c("Nonsense","Missense","Missense: Increased fitness\n(Obermeyer et al)")))%>%
  ggplot(aes(y=Gene)) +
  geom_vline(xintercept=1,linetype='dashed') +
  geom_pointrange(aes(x=Fold_change,xmin=minci,xmax=maxci,color=Mutation_type,fill=Mutation_type),position=position_dodge(width=0.4),size=1) +
  #facet_wrap(vars(Gene),nrow=3) +
  scale_x_log10(limits=c(0.00004,25),breaks=c(0.001,0.1,10),labels=c("0.001","0.1","10")) +
  theme_minimal() +
  xlab('Fold change in cluster size\nrelative to synonymous') +
  scale_fill_manual(name='Mutation type',values=c('#e71d36','#2EC4B6','grey30','lightgrey'))+
  scale_color_manual(name='Mutation type',values=c('#e71d36','#2EC4B6','grey30','lightgrey')) +
  ylab('') +
  theme(axis.text.y=element_text(size=11,color='black'))

ggsave("figs/fig3/clustSize_mutType_regression.pdf",dpi=300,width=6,height=4)






orf8_nb = glm.nb(n_descendants ~ mut_type + offset(log(time)),data=orf8,maxit=50)
orf8_noAlphaXBB = orf8[orf8$aa_mutations!='ORF8:Q27*' & orf8$aa_mutations!='ORF8:G8*',]
orf8_noAlphaXBB_nb = glm.nb(n_descendants ~ mut_type + offset(log(time)),data=orf8_noAlphaXBB,maxit=50)
orf8_obermeyer_nb = glm.nb(n_descendants ~ mut_type_obermeyer + offset(log(time)),data=orf8,maxit=100)
summary(orf8_nb)
summary(orf8_noAlphaXBB_nb)
summary(orf8_obermeyer_nb)
print(exp(cbind(Odds=coef(orf8_noAlphaXBB_nb),confint(orf8_noAlphaXBB_nb))))


print(exp(cbind(Odds=coef(orf8_obermeyer_nb),confint(orf8_obermeyer_nb))))
## Model fits better if pull out Obermeyer muts

orf1a_nb = glm.nb(n_descendants ~ mut_type + offset(log(time)),data=orf1a,maxit=100)
orf1a_obermeyer_nb = glm.nb(n_descendants ~ mut_type_obermeyer + offset(log(time)),data=orf1a,maxit=200)
summary(orf1a_nb)
summary(orf1a_obermeyer_nb)
print(exp(cbind(Odds=coef(orf1a_obermeyer_nb),confint(orf1a_obermeyer_nb))))
## obermeyer fits better, but I don't know what's up with the decreased ones


spike_nb = glm.nb(n_descendants ~ mut_type + offset(log(time)),data=spike,maxit=200)
spike_obermeyer_nb = glm.nb(n_descendants ~ mut_type_obermeyer + offset(log(time)),data=spike,maxit=200)
summary(spike_obermeyer_nb)
summary(spike_nb)
print(exp(cbind(Odds=coef(spike_nb),confint(spike_nb))))
## Yah spike is better without obermeyer








show_split = function(df) {
  counted = df %>% group_by(mut_type,fitness,mut_type_obermeyer) %>% count()
  return(counted)
}  

show_split(orf8)
show_split(orf1a)
show_split(spike)

library(jtools)
summ(orf8_pois,exp=T,confint = T)
summ(orf8_qpois,exp=T,confint=T)
summ(orf8_nb,exp=T,confint=T)
print(exp(cbind(Odds=coef(orf8_gamnb))))

spike_gamnb = mgcv::gam()


orf8_pois = run_pois(orf8)
1-pchisq(orf8_pois$deviance, orf8_pois$df.residual) 
## Unsurprisingly this model is way overdispersed

orf8_qpois = run_qpois(orf8)
1-pchisq(orf8_qpois$deviance, orf8_qpois$df.residual) 
## quasipoisson doesn't fix this problem

orf8_nb = run_nb(orf8)
orf8_nb = run_gamnb(orf8)

mgcv::gam.check(orf8_gamnb,type="pearson")

plot(ggeffects::ggpredict(orf8_gamnb), facets = TRUE)

orf8_ln = glm(leaf_count ~ time + time:mut_type, family=gaussian,data=orf8)
plot(residuals(orf8_ln) ~ fitted(orf8_ln,type="link"),xlab="fitted cluster size",ylab="Deviance residuals",pch=20,col="blue",log='y')
plot(orf8_ln$y ~ fitted(orf8_ln,type="link"),xlab="fitted cluster size",ylab="actual cluster size",pch=20,col="blue",log='y')
plot(residuals(orf8_ln) ~ orf8_ln$model$time,xlab="years in circulation",ylab="Deviance residuals",pch=20,col="blue",log='y')


summary(orf8_ln)
summ(orf8_ln)
qqnorm(ln.res)
qqline(ln.res)
plot(density(ln.res))


#spike_glm = run_glm(spike)
spike_glm2 = run_glm2(spike)

#orf1a_glm  = run_glm(orf1a)
orf1a_glm2 = run_glm2(orf1a)
