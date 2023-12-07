# Runs poisson regression on cluster size

# Load packages
library(tidyverse)
library(lubridate)
library(MASS)
library(VGAM)


#  Change working directory
setwd("~/Work/projects/covid/long-deletions")

# Read in data
load_data = function(path,ober,GENE, clades) {
  
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
    slice_max(rank,n=1) %>%
    left_join(clades, by =c("node_id" ="node")) %>%
    mutate(clade = if_else(clade == '21K (Omicron)'|clade == '21L (Omicron)'|clade == '21M (Omicron)','21K-M (Omicron)', clade)) %>%
    mutate(clade = if_else(clade == '21A (Delta)' | clade == '21I (Delta)'|clade=='21J (Delta)','21A,I,J (Delta)',clade))
  df$mut_type_obermeyer = as_factor(df$mut_type_obermeyer)
  df$mut_type_obermeyer = fct_relevel(df$mut_type_obermeyer,c('synonymous_no_increase'))
  df$mut_type = as_factor(df$mut_type)
  df$mut_type = fct_relevel(df$mut_type, c('synonymous','missense','nonsense'))
  
  return(df)
}

obermeyer = read_tsv('data/obermeyer_mutations.tsv')
clades = read_tsv('usher/trimmed/node_clades.tsv')
obermeyer$mutation = str_replace_all(obermeyer$mutation, "STOP", "*")
orf8 = load_data('usher/trimmed/clades_nested/ORF8_clades.tsv',obermeyer,'ORF8',clades)

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

coefficients =  function(fit,clade) {
  df = exp(cbind(Fold_change=coef(fit),confint(fit)))
  df <- cbind(Variable = rownames(df), df)
  df <- cbind(df,Clade=clade)
  return(as_tibble(df))
}

orf8_nb = nb_reg(orf8,50)
orf8_coef = coefficients(orf8_nb,'ORF8')

orf8_modified  = orf8 %>%
  mutate(effect = ifelse(mut_type=='synonymous','synonymous',ifelse(mut_type=='missense',)))
orf8_modified$variant = str_c(orf8_modified$clade,"_", orf8_modified$mut_type)

orf8_modified$variant = fct_relevel(orf8_modified$variant,'all_synonymous')
orf8_big = glm.nb(n_descendants ~ variant + offset(log(time)),data=orf8_modified,maxit=500)
summary(orf8_big)
full_results = coefficients(orf8_big,'joke')

full_results %>% 
  filter(Variable != '(Intercept)') %>%
  separate(Variable, c('clade','mutType'),sep='_') %>%
  mutate(clade = str_remove(clade,'variant')) %>%
  mutate(clade = if_else(clade %in% names(mapping), mapping[clade],clade)) %>%
  mutate(clade = factor(clade,levels=c('19A',
                                       '19B',
                                       '20A',
                                       '20B',
                                       '20C',
                                       '20D',
                                       '20E (EU1)',
                                       '20F',
                                       '20G',
                                       '20H (Beta,V2)',
                                       '20I (Alpha,V1)',
                                       '20J (Gamma,V3)',
                                       '21A (Delta)',
                                       '21A,I,J (Delta)',
                                       '21B (Kappa)',
                                       '21C (Epsilon)',
                                       '21D (Eta)',
                                       '21E (Theta)',
                                       '21F (Iota)',
                                       '21G (Lambda)',
                                       '21H (Mu)',
                                       '21I (Delta)',
                                       '21J (Delta)',
                                       '21K-M (Omicron)',
                                       '21K (Omicron)',
                                       '21L (Omicron)',
                                       '21M (Omicron)',
                                       '22A (Omicron)',
                                       "22A (BA.4)",
                                       '22B (Omicron)',
                                       "22B (BA.5)", 
                                       '22C (Omicron)',
                                       "22C (BA.2.12.1)",
                                       '22D (Omicron)',
                                       "22D (BA.2.75)", 
                                       '22E (Omicron)',
                                       "22E (BQ.1)",
                                       '22F (Omicron)',
                                       '22F (XBB)'))) %>%
  mutate(Fold_change = as.numeric(Fold_change)) %>%
  mutate(minci = as.numeric(`2.5 %`)) %>%
  mutate(maxci = as.numeric(`97.5 %`)) %>%
  mutate(clade = fct_rev(clade)) %>%
  ggplot(aes(y=clade)) +
  geom_vline(xintercept=1,linetype='dashed') +
  geom_pointrange(aes(x=Fold_change,xmin=minci,xmax=maxci,color=mutType,fill=mutType),position=position_dodge(width=0.4),size=1) +
  #facet_wrap(vars(Gene),nrow=3) +
  scale_x_log10(breaks=c(0.001,0.1,10),labels=c("0.001","0.1","10")) +
  theme_minimal() +
  coord_cartesian(xlim=c(0.0004,25))+
  xlab('Fold change in cluster growth rate\nrelative to ORF8 synonymous') +
  scale_fill_manual(name='Mutation type',values=c('#2EC4B6','#e71d36'),labels=c('Missense','Nonsense'))+
  scale_color_manual(name='Mutation type',values=c('#2EC4B6','#e71d36'),labels=c('Missense','Nonsense')) +
  #scale_fill_manual(name='Mutation type',values=c('#e71d36','#2EC4B6'))+
  #scale_color_manual(name='Mutation type',values=c('#e71d36','#2EC4B6')) +
  ylab('') +
  theme(axis.text.y=element_text(size=11,color='black'))

ggsave("figs/clustergrowthrate_ORF8_clade_background.pdf",dpi=300,width=6,height=8)
ggsave("figs/clustergrowthrate_ORF8_clade_background.jpg",dpi=300,width=6,height=8,bg='white')


orf8_nb_pos = nb_reg_pos(orf8,200)
unique(orf8$Clade)
for (clade in unique(orf8$clade)) {
  print(clade)
  if (clade != '20B') {
    filt = orf8[orf8$clade==clade,]
    reg = nb_reg(filt,500)
    coef = coefficients(reg,clade)
    
    if (exists("results")) {
      results = bind_rows(results,coef)
    } else {
      results = coef
    }
  }
}

format_results = function(df) {
  df$Fold_change = as.numeric(df$Fold_change)
  df$`2.5 %` = as.numeric(df$`2.5 %`)
  df$`97.5 %` = as.numeric(df$`97.5 %`)
  return(df)
}


results = format_results(results)

max(results[!is.na(results$`97.5 %`) & (results$Variable != '(Intercept)'),]$`97.5 %`)
mapping <- c("22F (Omicron)" = "22F (XBB)","22E (Omicron)"="22E (BQ.1)", "22D (Omicron)" = "22D (BA.2.75)", "22C (Omicron)" = "22C (BA.2.12.1)", "22B (Omicron)" = "22B (BA.5)", "22A (Omicron)" = "22A (BA.4)")
#%>%dplyr::select(clade) %>% distinct() %>% print(n=30) #
results %>%
  mutate(clade = if_else(Clade %in% names(mapping), mapping[Clade],Clade)) %>%
  filter(!grepl("Mu",clade)) %>%
  filter(!grepl("Lambda",clade)) %>%
  filter(!grepl("Theta",clade)) %>%
  filter(!grepl("Eta",clade)) %>%
  filter(!grepl("Kappa",clade)) %>%
  mutate(clade = factor(clade,levels=c('19A',
                                      '19B',
                                      '20A',
                                      '20B',
                                      '20C',
                                      '20D',
                                      '20E (EU1)',
                                      '20F',
                                      '20G',
                                      '20H (Beta,V2)',
                                      '20I (Alpha,V1)',
                                      '20J (Gamma,V3)',
                                      '21A (Delta)',
                                      '21A,I,J (Delta)',
                                      '21B (Kappa)',
                                      '21C (Epsilon)',
                                      '21D (Eta)',
                                      '21E (Theta)',
                                      '21F (Iota)',
                                      '21G (Lambda)',
                                      '21H (Mu)',
                                      '21I (Delta)',
                                      '21J (Delta)',
                                      '21K-M (Omicron)',
                                      '21K (Omicron)',
                                      '21L (Omicron)',
                                      '21M (Omicron)',
                                      '22A (Omicron)',
                                      "22A (BA.4)",
                                      '22B (Omicron)',
                                      "22B (BA.5)", 
                                      '22C (Omicron)',
                                      "22C (BA.2.12.1)",
                                      '22D (Omicron)',
                                      "22D (BA.2.75)", 
                                      '22E (Omicron)',
                                      "22E (BQ.1)",
                                      '22F (Omicron)',
                                      '22F (XBB)'))) %>%
  filter(Variable!='(Intercept)') %>%
  mutate(clade = fct_rev(clade)) %>%
  mutate(minci = replace_na(`2.5 %`,0.00004)) %>%
  mutate(maxci = replace_na(`97.5 %`,1600)) %>%
  ggplot(aes(y=clade)) +
  geom_vline(xintercept=1,linetype='dashed') +
  geom_pointrange(aes(x=Fold_change,xmin=minci,xmax=maxci,color=Variable,fill=Variable),position=position_dodge(width=0.4),size=1) +
  #facet_wrap(vars(Gene),nrow=3) +
  scale_x_log10(breaks=c(0.001,0.1,10),labels=c("0.001","0.1","10")) +
  theme_minimal() +
  #coord_cartesian(xlim=c(0.0004,25))+
  xlab('Fold change in cluster growth rate\nrelative to ORF8 synonymous') +
  scale_fill_manual(name='Mutation type',values=c('#2EC4B6','#e71d36'),labels=c('Missense','Nonsense'))+
  scale_color_manual(name='Mutation type',values=c('#2EC4B6','#e71d36'),labels=c('Missense','Nonsense')) +
  #scale_fill_manual(name='Mutation type',values=c('#e71d36','#2EC4B6'))+
  #scale_color_manual(name='Mutation type',values=c('#e71d36','#2EC4B6')) +
  ylab('') +
  theme(axis.text.y=element_text(size=11,color='black'))

ggsave("figs/clustergrowthrate_ORF8_clade.pdf",dpi=300,width=6,height=8)
ggsave("figs/clustergrowthrate_ORF8_clade.jpg",dpi=300,width=6,height=8,bg='white')

