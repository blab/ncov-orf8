library(tidyverse)
library(ggplot2)
library(ggridges)

setwd("~/Work/projects/covid/long-deletions/")

read_data = function(path) {
  df = read_tsv(path)
  df$mut_type = factor(df$mut_type,levels=c("synonymous","missense","nonsense","undoStop"))
  df <- df %>%
    mutate(daysPos = days_circulated + 1)
  return(df)
}


## Load data
orf8 = read_data('usher/trimmed/clades/ORF8_clades.tsv')
spike = read_data('usher/trimmed/clades/S_clades.tsv')
S1 = read_data('usher/trimmed/clades/S1_clades.tsv')
orf1a = read_data('usher/trimmed/clades/ORF1a_clades.tsv')

orf8_nested = read_data('usher/trimmed/clades_nested/ORF8_clades.tsv')
spike_nested = read_data('usher/trimmed/clades_nested/S_clades.tsv')
S1_nested = read_data('usher/trimmed/clades_nested/S1_clades.tsv')
orf1a_nested = read_data('usher/trimmed/clades_nested/ORF1a_clades.tsv')

obermeyer_muts = read_tsv('data/obermeyer_mutations.tsv')
obermeyer_muts$`R / R_A` = as.numeric(obermeyer_muts$`R / R_A`)

spike_merged = spike %>%
  left_join(obermeyer_muts, by = c("aa_mutations"="mutation")) %>%
  mutate(fitness = ifelse(`R / R_A 95% ci lower` > 1, 'fitness enhancing', ifelse(`R / R_A 95% ci upper` < 1, 'fitness decreasing', 'neutral')))

orf8 %>%
  group_by(mut_type) %>%
  add_count() %>%
  ungroup() %>%
  ggplot(aes(x=leaf_count,color=mut_type, y = (..count..)/n*100, group=mut_type)) +
  geom_bar() +
  scale_x_log10() + 
  scale_y_log10() 
  #geom_histogram()
  #geom_point(stat="bin", aes(y=..ndensity..), bins=20) +
  #geom_line(stat="bin", aes(y=..ndensity..), bins=20)
  #geom_boxplot() + 
  #stat_ecdf() +
  #geom_freqpoly() +
  ##geom_dotplot() +
