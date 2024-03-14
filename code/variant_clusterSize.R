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
  
  df = read_tsv(path)
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
    left_join(clades, by =c("node_id" ="node")) #%>%
    #mutate(clade = if_else(clade == '21K (Omicron)'|clade == '21L (Omicron)'|clade == '21M (Omicron)','21K-M (Omicron)', clade)) %>%
    #mutate(clade = if_else(clade == '21A (Delta)' | clade == '21I (Delta)'|clade=='21J (Delta)','21A,I,J (Delta)',clade))
  df$mut_type_obermeyer = as_factor(df$mut_type_obermeyer)
  df$mut_type_obermeyer = fct_relevel(df$mut_type_obermeyer,c('synonymous_no_increase'))
  df$mut_type = as_factor(df$mut_type)
  df$mut_type = fct_relevel(df$mut_type, c('synonymous','missense','nonsense', 'undoStop'))
  
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
  df <- cbind(df,Clade=clade,theta=fit$theta)
  return(as_tibble(df))
}

orf8_nb = nb_reg(orf8[orf8$leaf_count<100000,],50)
orf8_coef = coefficients(orf8_nb,'ORF8')
tree_theta = as.numeric(orf8_coef$theta[1])

clade_counts <- orf8 %>% group_by(clade,mut_type) %>% summarise(n=n())

clade_counts %>% 
  ggplot(aes(x=n,fill=mut_type)) +
  geom_histogram() +
  scale_x_log10() +
  facet_grid(.~ mut_type) +
  xlab('Number of mutations per clade')

## Same synonymous rate across entire tree
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

## Different synonymous rate by variant
for (clade in unique(orf8$clade)) {
  print(clade)
  if ((clade != '20B') & (clade != '21M (Omicron)') & (clade != '21E (Theta)') & (clade != '21G (Lambda)')) {
    filt = orf8[(orf8$clade==clade) &(orf8$mut_type!='undoStop'),]
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
write_tsv(results,'usher/trimmed/clusterGrowthRate_clade.tsv')

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

## What does it look like with the undo stop
filtAlpha = orf8[(orf8$clade=='20I (Alpha,V1)')&(orf8$leaf_count!=601725),]
fitAlpha = nb_reg(filtAlpha,maxit=200)
coefficients(fitAlpha,'Alpha')

filt22C = orf8[(orf8$clade=='22C (Omicron)')&(orf8$leaf_count!=8745),]
fit22C = nb_reg(filt22C,maxit=200)
coefficients(fit22C,'22C')

filtXBB = orf8[(orf8$clade=='22F (Omicron)')&(orf8$leaf_count!=116194),]
fitXBB = nb_reg(filtXBB,maxit=200)
coefficients(fitXBB,'XBB')

filtNoWinner = orf8[(orf8$clade != '20I (Alpha,V1)') & (orf8$clade != '22C (Omicron)') & (orf8$clade!='22F (Omicron)'),]
fitNoWinner = nb_reg(filtNoWinner,maxit=200)
coefficients(fitNoWinner,'without XBB, Alpha, Omicron')

filtNoExtr = orf8[(orf8$clade != '20I (Alpha,V1)') & (orf8$clade != '22C (Omicron)') & (orf8$clade!='22F (Omicron)') & (orf8$clade != '21J (Delta)') & (orf8$clade != '20E (EU1)'),]
fitNoExtr = nb_reg(filtNoExtr,maxit=200)
coefficients(fitNoExtr,'without XBB, Alpha, Omicron, 20E, and 21J')

## Exclude clades
for (clade in unique(orf8$clade)) {
  print(clade)
  filt = orf8[orf8$clade!=clade,]
  reg = nb_reg(filt,500)
  lab = paste('without', clade)
  coef = coefficients(reg,lab)
  
  if (exists("withoutResults")) {
    withoutResults = bind_rows(withoutResults,coef)
  } else {
    withoutResults = coef

  }
}

withoutResults <- format_results(withoutResults)

withoutResults %>%
  filter(Variable!='(Intercept)') %>%
  mutate(Clade = fct_rev(Clade)) %>%
  mutate(minci = replace_na(`2.5 %`,0.00004)) %>%
  mutate(maxci = replace_na(`97.5 %`,1600)) %>%
  ggplot(aes(y=Clade)) +
  geom_vline(xintercept=1,linetype='dashed') +
  geom_pointrange(aes(x=Fold_change,xmin=minci,xmax=maxci,color=Variable,fill=Variable),position=position_dodge(width=0.4),size=1) +
  #facet_wrap(vars(Gene),nrow=3) +
  scale_x_log10(breaks=c(0.001,0.1,10),labels=c("0.001","0.1","10")) +
  theme_minimal() +
  #coord_cartesian(xlim=c(0.0004,25))+
  xlab('Fold change in cluster growth rate\nrelative to ORF8 synonymous') +
  scale_fill_manual(name='Mutation type',values=c('#2EC4B6','#e71d36','#000000'),labels=c('Missense','Nonsense','UndoStop'))+
  scale_color_manual(name='Mutation type',values=c('#2EC4B6','#e71d36','#000000'),labels=c('Missense','Nonsense','UndoStop')) +
  #scale_fill_manual(name='Mutation type',values=c('#e71d36','#2EC4B6'))+
  #scale_color_manual(name='Mutation type',values=c('#e71d36','#2EC4B6')) +
  ylab('') +
  theme(axis.text.y=element_text(size=11,color='black'))
ggsave("figs/withoutClade_clusterGrowthRate_ORF8.pdf",dpi=300,width=6,height=8)

### Median cluster size
orf8 %>% group_by(mut_type,clade) %>% 
  filter(mut_type != 'undoStop') %>%
  #filter(mut_type != 'synonymous') %>%
  mutate(logSize = log10(leaf_count)) %>%
  summarise(median_size = median(leaf_count), per75 = quantile(leaf_count,c(0.75)), logmean = mean(logSize)) %>%
  ggplot(aes(x=clade,y=logmean,fill=mut_type,color=mut_type)) +
  #geom_point() +
  geom_bar(stat='identity',position=position_dodge()) +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_manual(values=c('#ff9f1c','#2EC4B6','#e71d36')) +
  scale_color_manual(values=c('#ff9f1c','#2EC4B6','#e71d36'))  + 
  #theme_minimal() +
  ylab('Mean of log10(clusterSize)')

orf8 %>% group_by(clade, mut_type) %>% summarise(n = n()) %>% filter(n<30) %>% filter(mut_type!='undoStop') %>% print(n=22)
