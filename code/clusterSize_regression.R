## Runs poisson regression on cluster size

# Load packages
library(tidyverse)
library(lubridate)
library(MASS)


#  Change working directory
setwd("~/Work/projects/covid/long-deletions")

# Read in data
load_data = function(path,ober,GENE) {
  df = read_tsv(path) %>% filter(mut_type!='undoStop')
  df$date_observed = ymd(df$date_observed)
  df$time = interval(df$date_observed, ymd(20230501)) %>% as.numeric("years")
  df <- df %>% 
    separate_longer_delim(c("aa_mutations","nt_mutations","codon_change"),delim=';') %>%
    left_join(ober, by=c("aa_mutations"="mutation")) %>%
    mutate(`Δ log R` = if_else(is.na(`Δ log R`),0,`Δ log R`)) %>%
    mutate(`Δ log R 95% ci upper` = if_else(is.na(`Δ log R 95% ci upper`),0,`Δ log R 95% ci upper`)) %>%
    mutate(`Δ log R 95% ci lower` = if_else(is.na(`Δ log R 95% ci lower`),0,`Δ log R 95% ci lower`)) %>%
    mutate(n_descendants = leaf_count - 1) %>%
    mutate(fitness = if_else(`Δ log R 95% ci lower`>0,'increased','no_increase')) %>%
    unite(mut_type_obermeyer,c("mut_type","fitness"),remove=FALSE) %>%
    filter(str_detect(aa_mutations,GENE)) %>%
    group_by(node_id,branch_muts_counts,gene,leaf_count,mut_count,days_circulated,date_observed,parent,cluster,time,n_descendants) %>%
    slice_max(rank,n=1)
  df$mut_type_obermeyer = as_factor(df$mut_type_obermeyer)
  df$mut_type_obermeyer = fct_relevel(df$mut_type_obermeyer,c('synonymous_no_increase'))
  df$mut_type = as_factor(df$mut_type)
  df$mut_type = fct_relevel(df$mut_type, c('synonymous','missense','nonsense'))
  return(df)
}

obermeyer = read_tsv('data/obermeyer_mutations.tsv')
obermeyer$mutation = str_replace_all(obermeyer$mutation, "STOP", "*")
orf8 = load_data('usher/trimmed/clades_nested/ORF8_clades.tsv',obermeyer,'ORF8')
orf1a = load_data('usher/trimmed/clades_nested/ORF1a_clades.tsv',obermeyer,'ORF1a')
spike = load_data('usher/trimmed/clades_nested/S_clades.tsv',obermeyer,'S')
orf1b = load_data('usher/trimmed/clades_nested/ORF1b_clades.tsv',obermeyer,'ORF1b')
orf3a = load_data('usher/trimmed/clades_nested/ORF3a_clades.tsv',obermeyer,'ORF3a')
e = load_data('usher/trimmed/clades_nested/E_clades.tsv',obermeyer,'E')
m = load_data('usher/trimmed/clades_nested/M_clades.tsv',obermeyer,'M')
orf6 = load_data('usher/trimmed/clades_nested/ORF6_clades.tsv',obermeyer,'ORF6')
orf7a = load_data('usher/trimmed/clades_nested/ORF7a_clades.tsv',obermeyer,'ORF7a')
orf7b = load_data('usher/trimmed/clades_nested/ORF7b_clades.tsv',obermeyer,'ORF7b')
n = load_data('usher/trimmed/clades_nested/N_clades.tsv',obermeyer,'N')
orf9b = load_data('usher/trimmed/clades_nested/ORF9b_clades.tsv',obermeyer,'ORF9b')


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



nb_reg = function(df,maxit,cutoff) {
  if(missing(cutoff)) {
    filt = df
  } else {
    filt = df[df$branch_muts_counts<cutoff,]
  }
  fit = glm.nb(n_descendants ~ mut_type + offset(log(time)),data=filt,maxit=maxit)
  print(summary(fit))
  return(fit)
}

nb_reg_obermeyer = function(df, maxit,cutoff) {
  if(missing(cutoff)) {
    filt = df
  } else {
    filt = df[df$branch_muts_counts<cutoff,]
  }
  fit = glm.nb(n_descendants ~ mut_type_obermeyer + offset(log(time)), data = filt, maxit = maxit)
  print(summary(fit))
  return(fit)
}

coefficients =  function(fit,gene){
  df = exp(cbind(Fold_change=coef(fit),confint(fit)))
  df <- cbind(Variable = rownames(df), df)
  df <- cbind(df,Gene=gene)
  return(as_tibble(df))
}

orf8_nb = nb_reg(orf8,50)

# Log likelihood 
pchisq(2 * (logLik(orf8_nb) - logLik(orf8_pois)), df = 4, lower.tail = FALSE)


orf8_coef = coefficients(orf8_nb,'ORF8')
orf1a_nb = nb_reg(orf1a,100)
orf1a_coef = coefficients(orf1a_nb,'ORF1a')
spike_nb = nb_reg(spike,200)
spike_coef = coefficients(spike_nb,'Spike')  
orf1b_nb = nb_reg(orf1b,100)
orf1b_coef = coefficients(orf1b_nb,'ORF1b')
orf3a_nb = nb_reg(orf3a,100)
orf3a_coef = coefficients(orf3a_nb,'ORF3a')  
e_nb = nb_reg(e[e$mut_type!='nonsense',],200)
e_coef = coefficients(e_nb,'E')  
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

results = bind_rows(orf1a_coef,orf1b_coef,spike_coef,orf3a_coef,e_coef,m_coef,orf6_coef,orf7a_coef,orf7b_coef,orf8_coef,n_coef,orf9b_coef)

orf1a_obermeyer_nb = nb_reg_obermeyer(orf1a, 100)
orf1a_obermeyer_coef = coefficients(orf1a_obermeyer_nb,'ORF1a')
orf1a_obermeyer_coef = cbind(orf1a_obermeyer_coef, Regression='withObermeyer')
spike_obermeyer_nb = nb_reg_obermeyer(spike, 200)
spike_obermeyer_coef = coefficients(spike_obermeyer_nb,'Spike')
spike_obermeyer_coef = cbind(spike_obermeyer_coef, Regression='withObermeyer')
orf8_obermeyer_nb = nb_reg_obermeyer(orf8,100)
orf8_obermeyer_coef = coefficients(orf8_obermeyer_nb,'ORF8')
orf8_obermeyer_coef = cbind(orf8_obermeyer_coef, Regression='withObermeyer')


results_small = bind_rows(orf1a_obermeyer_coef,spike_obermeyer_coef, orf8_obermeyer_coef,spike_coef,orf8_coef,orf1a_coef)
#results_small = bind_rows(orf1b_coef,spike_coef,orf8_coef,orf1a_coef,orf1b_obermeyer_coefwrit
write_delim(results_small,'results/clusterGrowthRate_small.tsv',delim='\t')

format_results = function(df) {
  df$Fold_change = as.numeric(df$Fold_change)
  df$`2.5 %` = as.numeric(df$`2.5 %`)
  df$`97.5 %` = as.numeric(df$`97.5 %`)
  df$Gene = as_factor(df$Gene)
  return(df)
}

results_small = format_results(results_small)

results_small$Gene = fct_relevel(results_small$Gene,c('ORF1b','Spike','ORF8'))
mapping = c("mut_type_obermeyermissense_increased"="Missense: Increased fitness\n(Obermeyer et al)", "mut_type_obermeyermissense_other"="Missense: Other",'mut_type_obermeyernonsense'="Nonsense","mut_typemissense"="Missense","mut_typenonsense"="Nonsense")


results = format_results(results)
## Palette
pal1 = c('#e71d36','#2EC4B6','#ff9f1c','#9E641B','#FFEBD2','#A2E1D8','#266962')
colorblindcheck::palette_check(pal1, plot = TRUE)


results_small %>%
  filter(Variable!='(Intercept)') %>%
  filter(Gene == 'Spike'|Gene =='ORF1a'|Gene =='ORF8') %>%
  #filter(Gene == 'Spike'|Gene =='ORF1b'|Gene =='ORF8') %>%
  mutate(minci = replace_na(`2.5 %`,0.001)) %>%
  mutate(maxci = replace_na(`97.5 %`,9.99)) %>%
  mutate(Mutation_type = mapping[Variable])%>% 
  mutate(Mutation_type = fct_relevel(Mutation_type, c("Nonsense","Missense","Missense: Increased fitness\n(Obermeyer et al)")))%>%
  #mutate(Mutation_type = fct_relevel(Mutation_type, c("Nonsense","Missense")))%>%
  ggplot(aes(y=Gene)) +
  geom_vline(xintercept=1,linetype='dashed') +
  geom_pointrange(aes(x=Fold_change,xmin=minci,xmax=maxci,color=Mutation_type,fill=Mutation_type),position=position_dodge(width=0.4),size=1) +
  #facet_wrap(vars(Gene),nrow=3) +
  scale_x_log10(limits=c(0.00004,25),breaks=c(0.001,0.1,10),labels=c("0.001","0.1","10")) +
  theme_minimal() +
  xlab('Fold change in cluster growth rate\nrelative to synonymous') +
  scale_fill_manual(name='Mutation type',values=c('#e71d36','#2EC4B6','grey30','lightgrey'))+
  scale_color_manual(name='Mutation type',values=c('#e71d36','#2EC4B6','grey30','lightgrey')) +
  #scale_fill_manual(name='Mutation type',values=c('#e71d36','#2EC4B6'))+
  #scale_color_manual(name='Mutation type',values=c('#e71d36','#2EC4B6')) +
  ylab('') +
  theme(axis.text.y=element_text(size=11,color='black'))

ggsave("figs/fig3/clustSize_mutType_regression.pdf",dpi=300,width=6,height=4)

max(results[!is.na(results$`97.5 %`) & (results$Variable != '(Intercept)'),]$`97.5 %`)

results %>%
  filter(Variable!='(Intercept)') %>%
  filter(Gene!='ORF9b') %>%
  mutate(minci = replace_na(`2.5 %`,0.00004)) %>%
  mutate(maxci = replace_na(`97.5 %`,100)) %>%
  mutate(Mutation_type = mapping[Variable])%>% 
  mutate(Mutation_type = fct_relevel(Mutation_type, c("Nonsense","Missense")))%>%
  mutate(Gene = fct_rev(Gene)) %>%
  ggplot(aes(y=Gene)) +
  geom_vline(xintercept=1,linetype='dashed') +
  geom_pointrange(aes(x=Fold_change,xmin=minci,xmax=maxci,color=Mutation_type,fill=Mutation_type),position=position_dodge(width=0.4),size=1) +
  #facet_wrap(vars(Gene),nrow=3) +
  scale_x_log10(breaks=c(0.001,0.1,10),labels=c("0.001","0.1","10")) +
  theme_minimal() +
  coord_cartesian(xlim=c(0.0004,25))+
  xlab('Fold change in cluster growth rate\nrelative to synonymous') +
  scale_fill_manual(name='Mutation type',values=c('#e71d36','#2EC4B6'))+
  scale_color_manual(name='Mutation type',values=c('#e71d36','#2EC4B6')) +
  #scale_fill_manual(name='Mutation type',values=c('#e71d36','#2EC4B6'))+
  #scale_color_manual(name='Mutation type',values=c('#e71d36','#2EC4B6')) +
  ylab('') +
  theme(axis.text.y=element_text(size=11,color='black'))

ggsave("figs/supplemental/clustergrowthrate_gene.pdf",dpi=300,width=6,height=6)
ggsave("figs/supplemental/clustergrowthrate_gene.jpg",dpi=300,width=6,height=6,bg='white')










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



orf8_nb = nb_reg(orf8,50,2)
orf8_coef = coefficients(orf8_nb,'ORF8')
orf1a_nb = nb_reg(orf1a,100,2)
orf1a_coef = coefficients(orf1a_nb,'ORF1a')
spike_nb = nb_reg(spike,200,2)
spike_coef = coefficients(spike_nb,'Spike')  
orf1b_nb = nb_reg(orf1b,100,2)
orf1b_coef = coefficients(orf1b_nb,'ORF1b')
orf3a_nb = nb_reg(orf3a,100,2)
orf3a_coef = coefficients(orf3a_nb,'ORF3a')  
e_nb = nb_reg(e[e$mut_type!='nonsense',],200,2)
e_coef = coefficients(e_nb,'E')  
m_nb = nb_reg(m,300,2)
m_coef = coefficients(m_nb,'M')  
orf6_nb = nb_reg(orf6,300,2)
orf6_coef = coefficients(orf6_nb,'ORF6')
orf7a_nb = nb_reg(orf7a,100,2)
orf7a_coef = coefficients(orf7a_nb,'ORF7a')
orf7b_nb = nb_reg(orf7b,300,2)
orf7b_coef = coefficients(orf7b_nb,'ORF7b')  
n_nb = nb_reg(n,200,2)
n_coef = coefficients(n_nb,'N')  
orf9b_nb = nb_reg(orf9b,200,2)
orf9b_coef = coefficients(orf9b_nb,'ORF9b')  

results = bind_rows(orf1a_coef,orf1b_coef,spike_coef,orf3a_coef,e_coef,m_coef,orf6_coef,orf7a_coef,orf7b_coef,orf8_coef,n_coef,orf9b_coef)

orf1a_obermeyer_nb = nb_reg_obermeyer(orf1a, 100)
orf1a_obermeyer_coef = coefficients(orf1a_obermeyer_nb,'ORF1a')
orf1b_obermeyer_nb = nb_reg_obermeyer(orf1b, 100)
orf1b_obermeyer_coef = coefficients(orf1b_obermeyer_nb,'ORF1b')

results_small = bind_rows(orf1a_obermeyer_ceof,spike_coef,orf8_coef,orf1a_coef)
#results_small = bind_rows(orf1b_coef,spike_coef,orf8_coef,orf1a_coef,orf1b_obermeyer_coef)


results_small = format_results(results_small)

results_small$Gene = fct_relevel(results_small$Gene,c('ORF1a','Spike','ORF8'))


results = format_results(results)
## Palette
pal1 = c('#e71d36','#2EC4B6','#ff9f1c','#9E641B','#FFEBD2','#A2E1D8','#266962')
colorblindcheck::palette_check(pal1, plot = TRUE)


results_small %>%
  filter(Variable!='(Intercept)') %>%
  filter(Gene == 'Spike'|Gene =='ORF1a'|Gene =='ORF8') %>%
  #filter(Gene == 'Spike'|Gene =='ORF1b'|Gene =='ORF8') %>%
  mutate(minci = replace_na(`2.5 %`,0.001)) %>%
  mutate(maxci = replace_na(`97.5 %`,9.99)) %>%
  mutate(Mutation_type = mapping[Variable])%>% 
  mutate(Mutation_type = fct_relevel(Mutation_type, c("Nonsense","Missense","Missense: Increased fitness\n(Obermeyer et al)")))%>%
  #mutate(Mutation_type = fct_relevel(Mutation_type, c("Nonsense","Missense")))%>%
  ggplot(aes(y=Gene)) +
  geom_vline(xintercept=1,linetype='dashed') +
  geom_pointrange(aes(x=Fold_change,xmin=minci,xmax=maxci,color=Mutation_type,fill=Mutation_type),position=position_dodge(width=0.4),size=1) +
  #facet_wrap(vars(Gene),nrow=3) +
  scale_x_log10(limits=c(0.00004,25),breaks=c(0.001,0.1,10),labels=c("0.001","0.1","10")) +
  theme_minimal() +
  xlab('Fold change in cluster growth rate\nrelative to synonymous') +
  scale_fill_manual(name='Mutation type',values=c('#e71d36','#2EC4B6','grey30','lightgrey'))+
  scale_color_manual(name='Mutation type',values=c('#e71d36','#2EC4B6','grey30','lightgrey')) +
  #scale_fill_manual(name='Mutation type',values=c('#e71d36','#2EC4B6'))+
  #scale_color_manual(name='Mutation type',values=c('#e71d36','#2EC4B6')) +
  ylab('') +
  theme(axis.text.y=element_text(size=11,color='black'))

#ggsave("figs/fig3/clustSize_mutType_regression.pdf",dpi=300,width=6,height=4)

results %>%
  filter(Variable!='(Intercept)') %>%
  filter(Gene!='ORF9b') %>%
  mutate(minci = replace_na(`2.5 %`,0.00004)) %>%
  mutate(maxci = replace_na(`97.5 %`,100)) %>%
  mutate(Mutation_type = mapping[Variable])%>% 
  mutate(Mutation_type = fct_relevel(Mutation_type, c("Nonsense","Missense")))%>%
  mutate(Gene = fct_rev(Gene)) %>%
  ggplot(aes(y=Gene)) +
  geom_vline(xintercept=1,linetype='dashed') +
  geom_pointrange(aes(x=Fold_change,xmin=minci,xmax=maxci,color=Mutation_type,fill=Mutation_type),position=position_dodge(width=0.4),size=1) +
  #facet_wrap(vars(Gene),nrow=3) +
  scale_x_log10(breaks=c(0.001,0.1,10),labels=c("0.001","0.1","10")) +
  theme_minimal() +
  coord_cartesian(xlim=c(0.0004,25))+
  xlab('Fold change in cluster growth rate\nrelative to synonymous') +
  scale_fill_manual(name='Mutation type',values=c('#e71d36','#2EC4B6'))+
  scale_color_manual(name='Mutation type',values=c('#e71d36','#2EC4B6')) +
  #scale_fill_manual(name='Mutation type',values=c('#e71d36','#2EC4B6'))+
  #scale_color_manual(name='Mutation type',values=c('#e71d36','#2EC4B6')) +
  ylab('') +
  theme(axis.text.y=element_text(size=11,color='black'))

#ggsave("figs/supplemental/clustergrowthrate_gene.pdf",dpi=300,width=6,height=6)
#ggsave("figs/supplemental/clustergrowthrate_gene.jpg",dpi=300,width=6,height=6,bg='white')















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
