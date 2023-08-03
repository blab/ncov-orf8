library(tidyverse)
library(edgee)
setwd("~/Work/projects/covid/long-deletions/")

## Read in data
df = read_tsv('usher/trimmed/orf8clades.tsv')
noSingleton <- df %>%
  filter(leaf_count>1)

## Code for running Edgeworth
runEdgeworth <- function(df, type1, type2) {
  filt = df %>%
    filter(mut_type == type1 | mut_type == type2) %>%
    mutate(mut_type = ifelse(mut_type == 'synonymous',0,ifelse(mut_type=='missense', 1, 2))) %>%
    select(mut_type,leaf_count,days_circulated)
  
  mat = t(as.matrix(filt))
  result = empEdge(mat[2:3,],a=mat[1,],side='right',type='Welch')
  return(result)
}

## Missense vs. synonymous
runEdgeworth(df,'missense','synonymous')
runEdgeworth(noSingleton,'missense','synonymous')

## Nonsense vs. synonymous
runEdgeworth(df,'nonsense','synonymous')
runEdgeworth(noSingleton,'nonsense','synonymous')

## Nonsense vs. missense
runEdgeworth(df,'nonsense','missense')
runEdgeworth(noSingleton,'nonsense','missense')
