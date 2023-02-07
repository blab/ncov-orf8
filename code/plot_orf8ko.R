library(tidyverse)
library(lubridate)

# Read in SARS2 ORF8 dataset
ko <- read_tsv("washington.orf8ko.tsv")
metadata <- read_tsv("data/2022-07-26.gisaid.washington.metadata.tsv")

freq_gaps <- count(ko, gap)

## Plot distributon of gaps
ko %>%
  ggplot(aes(gap)) +
  geom_density() +
  theme_minimal() +
  xlab("Number of gaps in Orf8") +
  ylab("Density in Washington SARS2 sequences") +
  annotate(geom="text", x=0, y=0.63, label = "0 gaps") +
  annotate(geom="text", x = 16, y= 0.25, label = "6 gaps")

ggsave("Density_gaps.pdf")
  
ko %>%
  filter(gap > 6) %>%
  ggplot(aes(gap)) + 
  geom_histogram(binwidth=5) + theme_minimal() +
  xlab("Number of gaps in Orf 8") + 
  ylab("Count Washington SARS2 sequences") +
  scale_x_continuous(expand= c(0,1), limits = c(6, NA)) +
  annotate(geom="text", x=200, y= 125, label = "Only samples with more than 6 gaps", size=6)

ggsave("hist_7gaps.pdf")

## Plot distribution of N's
ko %>%
  ggplot(aes(Ns)) +
  geom_histogram() +
  theme_minimal() +
  xlab("Number of N's in Orf 8") +
  ylab("Count in Washington SARS2 sequences")

ggsave("hist_ns.pdf")

freq_ns <- count(ko,Ns)

## Plot distribution of protein lengths
ko %>%
  ggplot(aes(protein_length)) +
  geom_histogram() +
  theme_minimal() +
  xlab("Length of translated Orf 8") +
  ylab("Count in Washington samples")

ggsave("hist_length.pdf")

freq_lengths <- count(ko, protein_length)

## Call everything with more than 6 gaps, a protein length of <= 100 or more than 100 Ns an Orf 8 ko

 ko <- ko %>%
  mutate(orf8ko = ifelse(gap > 6 | protein_length <= 100 | Ns > 100, "yes", "no"))

 count(ko, orf8ko)

df <- metadata %>%
  select(strain, date, division, location, Nextstrain_clade, pango_lineage) %>%
  right_join(ko, by = "strain")

## Adds month and date
df <- df %>%
  mutate(month = month(date), year = year(date))

## What is distribution of ORF8 knockouts by time?

df %>%
  ggplot(aes(x=date, fill =orf8ko, color=orf8ko)) +
  geom_histogram(position="stack",bins=50) +
  theme_minimal() +
  xlab("Date") +
  ylab("Washington SARS2 sequences") +
  scale_fill_discrete(name="Orf8 KO", labels = c("No", "Yes")) +
  scale_color_discrete(name="Orf8 KO", labels = c("No", "Yes"))

ggsave("Orf8ko_time.pdf")


df %>%
  ggplot(aes(y=Nextstrain_clade,fill=orf8ko)) +
  geom_bar() +
  theme_minimal() +
  scale_fill_discrete(name="Orf8 KO", labels = c("No", "Yes")) +
  #scale_color_discrete(name="Orf8 KO", labels = c("No", "Yes")) +
  ylab("Nextstrain clade") +
  xlab("Count")

ggsave("Of8ko_clade.pdf")

alpha_ko <- df %>%
  filter(orf8ko=="yes") %>%
  filter(Nextstrain_clade=="20I (Alpha, V1)")


alpha_ko %>% 
  count(protein_length)

omicron_ko <- df %>%
  filter(orf8ko=="yes") %>%
  filter(str_detect(Nextstrain_clade, "Omicron"))

omicron_ko %>%
  count(protein_length<100)

## Get time window of ORF8 KO
df %>%
  filter(!is.na(date)) %>%
  group_by(orf8ko) %>%
  summarise(max_date = max(date), min_date = min(date))

count_time_orf8ko <- df %>%
  group_by(month, year, orf8ko) %>%
  count() %>%
  arrange(year)
