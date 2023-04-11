library(tidyverse)
library(lubridate)

setwd("~/Work/projects/covid/long-deletions")

# Read in SARS2 ORF8 dataset
df <- read_tsv("results/gisaid.washington.june20-july22.orf8ko.meta.tsv")
clust <-read_tsv("clusters/geneClusters/WA_20K_ORF8_clusters.tsv")

## What is distribution of ORF8 knockouts by time?

df %>%
  filter(!is.na(ORF8_ko)) %>%
  ggplot(aes(x=date, fill =ORF8_ko, color=ORF8_ko)) +
  geom_histogram(position="stack",bins=50) +
  theme_minimal() +
  xlab("Date") +
  ylab("Washington\nSARS-CoV-2\nsequences") +
  scale_fill_manual(name="Orf8 KO", labels = c("No", "Yes"),values = c('#2a9d8f', '#ce4257')) +
  scale_color_manual(name="Orf8 KO", labels = c("No", "Yes"),values = c('#2a9d8f', '#ce4257')) +
  ylim(0,7000)

ggsave("figs/Orf8ko_time_updated.pdf",dpi=300,width=6.5,height=2)


df %>%
  filter(!is.na(ORF8_ko)) %>%
  ggplot(aes(y=Nextstrain_clade, fill =ORF8_ko, color=ORF8_ko)) +
  geom_histogram(stat="count",position='fill',bins=50) +
  theme_minimal() +
  xlab("Proportion") +
  ylab("") +
  scale_fill_manual(name="Orf8 KO", labels = c("No", "Yes"),values = c('#2a9d8f', '#ce4257')) +
  scale_color_manual(name="Orf8 KO", labels = c("No", "Yes"),values = c('#2a9d8f', '#ce4257')) +
  theme(legend.position="none")

ggsave("figs/Of8ko_clade_updated.pdf",dpi=300,width=3,height=4.5)


clust %>%
  group_by(cluster,koType) %>%
  count() %>%
  filter(n>1) %>%
  ggplot(aes(x=koType,y=n)) +
  geom_violin(fill='#ce4257') +
  #geom_dotplot(binaxis='y',stackdir='center') +
  scale_y_log10() +
  ylab('Cluster size') +
  xlab('KO type') +
  scale_x_discrete(labels=c('Large deletion','Early stop codon')) +
  theme_minimal()

ggsave("figs/Of8ko_clustSize.pdf",dpi=300,width=3,height=2)

