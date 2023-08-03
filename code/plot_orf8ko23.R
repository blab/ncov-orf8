library(tidyverse)
library(lubridate)
library(colorblindcheck)

setwd("~/Work/projects/covid/long-deletions")

# Read in SARS2 ORF8 dataset
df <- read_tsv("results/gisaid.washington_ko_meta.tsv")
clust <-read_tsv("nextstrain_helper/results/WA_20k/clusters/clusters_ORF8.tsv")
df <- df %>%
  mutate(ORF8_koType = ifelse(is.na(ORF8_koType),'None', ORF8_koType))

df$ORF8_koType = as_factor(df$ORF8_koType)
df$ORF8_koType = fct_relevel(df$ORF8_koType, c("None","earlyStop","BigDeletion"))


# What % are KO?
nrow(df[df$ORF8_ko=="Yes" & df$coverage>=0.95,])

nrow(df[df$ORF8_ko=="Yes" & df$coverage>=0.95,])/nrow(df[df$coverage>=0.95,])

length(df)

## Palette
pal1 = c('#29335C','#C9B6BE','#707070')
pal2 = c('#3D4C89','#80A26B','#C6C6C6')
palette_check(pal1, plot = TRUE)
palette_check(pal2, plot = TRUE)

## What is distribution of ORF8 knockouts by time?


df %>%
  filter(coverage>=0.95) %>%
  ggplot(aes(x=date, fill =ORF8_koType, color=ORF8_koType)) +
  geom_histogram(position="stack",binwidth=7) +
  theme_minimal() +
  xlab("Date") +
  ylab("Washington State \nSARS-CoV-2 sequences") +
  scale_fill_manual(name="ORF8 KO", labels = c("None", "Premature Stop", "Large deletion"),values = c('#C6C6C6','#3D4C89','#80A26B')) +
  scale_color_manual(name="ORF8 KO", labels = c("None", "Premature Stop", "Large deletion"),values = c('#C6C6C6','#3D4C89','#80A26B')) +
  ylim(0,4000) +
  theme(legend.position=c(0.2,0.75))

ggsave("figs/fig1/Orf8ko_time.pdf",dpi=300,width=5,height=2.5)


## What is the distribution of ORF8KO by clades
clades = c("23B (Omicron)" = "23B (XBB.1.16)", "23A (Omicron)" = "23A (XBB.1.5)", "22F (Omicron)" = "22F (XBB)","22E (Omicron)"="22E (BQ.1)", "22D (Omicron)" = "22D (BA.2.75)", "22C (Omicron)" = "22C (BA.2.12.1)", "22B (Omicron)" = "22B (BA.5)", "22A (Omicron)" = "22A (BA.4)")

no_alphaXBB <- df %>% 
  filter(coverage>=0.95) %>% 
  mutate(clade = if_else(Nextstrain_clade %in% names(clades), clades[Nextstrain_clade], Nextstrain_clade)) %>%
  filter(!str_detect(clade,'XBB')) %>%
  filter(!str_detect(clade, 'Alpha'))

nrow(no_alphaXBB[no_alphaXBB$ORF8_ko=="Yes",])/nrow(no_alphaXBB)

df %>%
  filter(coverage>=0.95) %>%
  mutate(clade = if_else(Nextstrain_clade %in% names(clades), clades[Nextstrain_clade], Nextstrain_clade)) %>%  
  ggplot(aes(y=clade, fill =ORF8_koType, color=ORF8_koType)) +
  geom_histogram(stat="count",position='fill',bins=50) +
  theme_minimal() +
  xlab("Proportion") +
  ylab("") +
  scale_fill_manual(name="ORF8 KO", labels = c("None", "Premature Stop", "Large deletion"),values = c('#C6C6C6','#3D4C89','#80A26B')) +
  scale_color_manual(name="ORF8 KO", labels = c("None", "Premature Stop", "Large deletion"),values = c('#C6C6C6','#3D4C89','#80A26B')) +
  theme(legend.position="none")

ggsave("figs/fig1/Orf8ko_clade.pdf",dpi=300,width=2.8,height=4)

## What is the distribution of KO by type?
nrow(df[df$ORF8_ko=="Yes" & df$coverage>=0.95,])

nrow(df[df$ORF8_ko=="Yes" & df$coverage>=0.95 & df$ORF8_koType=='BigDeletion',])/nrow(df[df$ORF8_ko=="Yes" & df$coverage>=0.95,])

clust %>%
  group_by(cluster,ORF8_koType) %>%
  count() %>%
  filter(n>1) %>%
  ggplot(aes(x=ORF8_koType,y=n,color=ORF8_koType,fill=ORF8_koType)) +
  geom_violin() +
  scale_fill_manual(name="ORF8 KO", labels = c("Premature Stop", "Large deletion"),values = c('#80A26B','#3D4C89')) +
  scale_color_manual(name="ORF8 KO", labels = c("Premature Stop", "Large deletion"),values = c('#80A26B','#3D4C89')) +
  #geom_dotplot(binaxis='y',stackdir='center') +
  scale_y_log10() +
  ylab('Cluster size') +
  xlab('') +
  scale_x_discrete(labels=c('Large\ndeletion','Premature\nstop')) +
  theme_minimal() +
  theme(legend.position="none") +
  theme(axis.text.x=element_text(size=11,color='black'))

ggsave("figs/fig1/Orf8ko_clustSize.pdf",dpi=300,width=2.5,height=2)



clust <- clust %>%
  separate(ORF8_deletions, c("start","end"),sep=':',remove=FALSE) %>%
  mutate(start = ifelse(start == "None",NA,start)) %>%
  group_by(cluster) %>%
  add_count() %>%
  ungroup() %>%
  separate(ORF8_misStops,c("bp","codon"),sep=':',remove=FALSE) %>%
  separate(bp,c("start_bp","end_bp"),sep='-',remove=TRUE)
clust$start_bp = as.numeric(clust$start_bp)
clust$protein_size = (clust$start_bp + 1 - 27894)/3
clust$start = as.numeric(clust$start)
clust$end = as.numeric(clust$end)
clust$deletionSize = clust$end-clust$start

clust.grouped <- clust %>% group_by(cluster,ORF8_koType,deletionSize,protein_size) %>% count()
clust.grouped %>%
  filter(ORF8_koType=='bigDeletion') 

clust.grouped %>%
  filter(ORF8_koType=='earlyStop')

clust.grouped %>%
  filter(n>1)

mean(clust.grouped[clust.grouped$ORF8_koType=='bigDeletion',]$n)
mean(clust.grouped[clust.grouped$ORF8_koType=='earlyStop',]$n)

wilcox.test(clust.grouped[clust.grouped$ORF8_koType=='bigDeletion',]$n,clust.grouped[clust.grouped$ORF8_koType=='earlyStop',]$n)


cor.test(clust.grouped[(clust.grouped$ORF8_koType=='bigDeletion')&(clust.grouped$n>1),]$deletionSize,clust.grouped[(clust.grouped$ORF8_koType=='bigDeletion')&(clust.grouped$n>1),]$n)
cor.test(clust.grouped[(clust.grouped$ORF8_koType=='bigDeletion'),]$deletionSize,clust.grouped[(clust.grouped$ORF8_koType=='bigDeletion'),]$n)

dsize_clust <- clust.grouped %>%
  filter(ORF8_koType=='bigDeletion') %>%
  filter(n>1) %>%
  ggplot(aes(x=deletionSize,y=n)) +
  geom_point() +
  theme_minimal() +
  xlab('Deletion size') +
  ylab('Size of cluster') +
  geom_smooth(method="lm") +
  annotate("text",x=100,y=13,label="Pearson's: 0.47\np-value=0.012",size=3)+
  theme(text=element_text(family='Helvetica'))
dsize_clust
ggsave("figs/supplemental/deletionSize_clustSize.pdf",dpi=300,width=3,height=2)

cor.test(clust.grouped[(clust.grouped$ORF8_koType=='earlyStop')&(clust.grouped$n>1),]$protein_size,clust.grouped[(clust.grouped$ORF8_koType=='earlyStop')&(clust.grouped$n>1),]$n)
cor.test(clust.grouped[(clust.grouped$ORF8_koType=='earlyStop'),]$protein_size,clust.grouped[(clust.grouped$ORF8_koType=='earlyStop'),]$n)

psize_clust <- clust.grouped %>%
  filter(ORF8_koType=='earlyStop') %>%
  filter(n>1) %>%
  ggplot(aes(x=protein_size,y=n)) +
  geom_point() +
  theme_minimal() +
  xlab('Truncated protein size') +
  ylab('Size of cluster') +
  geom_smooth(method="lm") +
  scale_y_log10() + 
  annotate("text",x=90,y=500,label="Pearson's: -0.14\np-value=0.49",size=3) +
  theme(text=element_text(family='Helvetica'))
psize_clust
ggsave("figs/supplemental/proteinSize_clustSize.pdf",dpi=300,width=3,height=2)

biggestProtein = max(clust.grouped[!is.na(clust.grouped$protein_size) & clust.grouped$n>1,]$protein_size)

clust.grouped %>%
  filter(n>1) %>%
  filter(ORF8_koType=='earlyStop') %>%
  ggplot(aes(x=protein_size)) +
  geom_histogram(bins=biggestProtein) +
  theme_minimal() +
  ylab('Number of clusters') +
  xlab('Truncated protein size')

biggestDeletion = max(clust.grouped[!is.na(clust.grouped$deletionSize) & clust.grouped$n>1,]$deletionSize)

clust.grouped %>%
  filter(n>1) %>%
  filter(ORF8_koType=='bigDeletion') %>%
  ggplot(aes(x=deletionSize)) +
  geom_histogram(bins=biggestDeletion) +
  theme_minimal() +
  ylab('Number of clusters') +
  xlab('Deletion size')

deletion_dist <- clust.grouped %>%
  filter(n>1) %>%
  filter(ORF8_koType=='bigDeletion') %>%
  ggplot(aes(x=deletionSize)) +
  geom_boxplot(fill='#80A26B') +
  xlab('Deletion size') +
  theme_minimal()+
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank(), panel.grid.major = element_blank(),
         text=element_text(family='Helvetica'))

protein_dist <- clust.grouped %>%
  filter(n>1) %>%
  filter(ORF8_koType=='earlyStop') %>%
  ggplot(aes(x=protein_size)) +
  geom_boxplot(fill='#3D4C89') +
  theme_minimal()+
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank(), panel.grid.major = element_blank(),text=element_text(family='Helvetica')) +
  xlab('Truncated protein size')


library(ggpubr)
figure <- ggarrange(dsize_clust,psize_clust, deletion_dist, protein_dist,
                    labels = c("A", "B", "C","D"),
                    ncol = 2, nrow = 2,heights=c(3,1))

figure
ggsave("figs/supplemental/proteinSize_deletionSize_clusterSize.pdf",dpi=300,width=6,height=3)
ggsave("figs/supplemental/proteinSize_deletionSize_clusterSize.jpg",dpi=300,width=6,height=3,bg='white')


nrow(clust.grouped[clust.grouped$n>1 & clust.grouped$ORF8_koType=='earlyStop' & clust.grouped$protein_size<= 26,])/nrow(clust.grouped[clust.grouped$n>1 & clust.grouped$ORF8_koType=='earlyStop',])

clust.wa <- clust %>% inner_join(df,by=c("strain"))

adj <- clust.wa %>% filter(n>1) %>%
  filter(ORF8_koType.x=='bigDeletion') %>%
  group_by(ORF7a_ko,ORF7b_ko,N_ko,deletionSize,cluster)%>%
  count() %>%
  arrange(cluster)

clust.grouped %>%
  filter(deletionSize>300) %>%
  filter(n>1)

