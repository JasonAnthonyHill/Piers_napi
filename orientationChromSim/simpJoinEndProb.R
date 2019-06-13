## Simplified end joining probability
library(dplyr)
library(data.table)
library(ggplot2)

source('syntRuns.R')

# Record highest number blocks (end blocks) on each genome
BmPnSyntRuns <- BmoriPnapiSynt %>% 
  arrange(BmChr, BmGeneStart) %>% 
  group_by(BmChr) %>% 
  mutate(lastBmBlk = max(BmRunID)) %>%
  ungroup() %>% 
  arrange(PnChr, PnGeneStart) %>% 
  group_by(PnChr) %>% 
  mutate(lastPnBlk = max(runid)) %>% 
  ungroup()

# Assign direction to synt blocks on the Pn chromosomes
blockdir <- BmPnSyntRuns %>% 
  arrange(PnChr, PnGeneStart) %>% 
  group_by(sbid) %>% 
  mutate(BmBlkDir = (lead(BmGeneStart) - BmGeneStart)) %>%
  summarise(sbum = sum(BmBlkDir, na.rm = T)) %>% 
  mutate(sbForward = ifelse(sbum <= 0, 0, 1)) %>% 
  select(sbid, sbForward)

slimSyntRuns <- BmPnSyntRuns %>%
  select(BmChr, PnChr, sbid, runid, BmRunID, lastBmBlk, lastPnBlk) %>% 
  group_by(sbid) %>% 
  slice(1) %>% 
  left_join(blockdir, by = 'sbid')

## Call blocks by their location in Bm and Pn genomes
slimSyntRuns <- slimSyntRuns %>%
  rename(PnRunID = runid) %>% 
  mutate(BmPos = ifelse(BmRunID == 1, 'b',
                        ifelse(BmRunID == lastBmBlk, 'e', 'i'))) %>% 
  mutate(PnPos = ifelse(PnRunID == 1, 'b',
                        ifelse(PnRunID == lastPnBlk, 'e', 'i')))

## Encode end of chromosome matching
slimSyntRuns <- slimSyntRuns %>% 
  mutate(matchEnd = ifelse(BmPos == 'i' | PnPos == 'i', 0,
                          ifelse(PnPos == 'b' & BmPos == 'b' & sbForward == 1, 1,
                                 ifelse(PnPos == 'e' & BmPos == 'e' & sbForward == 1, 1,
                                        ifelse(PnPos == 'b' & BmPos == 'e' & sbForward == 0, 1,
                                               ifelse(PnPos == 'e' & BmPos == 'b' & sbForward == 0, 1,0))))))

ObservedEnds <- sum(slimSyntRuns$matchEnd)

## Calculate chromosome end maintainence for random genome joins
## Function to rearrange syneny data frame by random draw
reordChrom2 <- function(df.in){
  df.in <- BmoriPnapiSynt
  sbids <- df.in$sbid %>% unique()
  ord <- sample(sbids, replace = FALSE, length(sbids))
  gn.df <- split(ord, ceiling(seq_along(ord)/4)) %>% 
    plyr::ldply(data.frame, .id = 'FuseChr') %>% 
    rename(sbid = X..i..) %>% 
    mutate(orient = sample(c('f','r'), length(sbids), replace = T)) %>% 
    right_join(df.in, by = 'sbid') %>% 
    mutate(sbid = factor(sbid, levels = ord)) %>% 
    group_by(sbid) %>%  
    arrange(sbid, ifelse(orient == "f", BmGeneStart, -BmGeneStart)) %>%
    ungroup(sbid) %>%
    group_by(FuseChr) %>% 
    mutate(FuseRunID = rleid(BmChr))
  return(gn.df)
}

countEnds <- function(df2.in){
  BmFuseSyntRuns <- df2.in %>% 
    arrange(BmChr, BmGeneStart) %>% 
    group_by(BmChr) %>% 
    mutate(lastBmBlk = max(BmRunID)) %>%
    ungroup() %>% 
    group_by(FuseChr) %>% 
    mutate(lastFuseBlk = max(FuseRunID)) %>% 
    ungroup()
  
  # Assign direction to synt blocks on the Fuse chromosomes
  blockdir <- BmFuseSyntRuns %>% 
    arrange(FuseChr, ifelse(orient == 'f', PnGeneStart, -PnGeneStart)) %>% 
    group_by(sbid) %>% 
    mutate(BmBlkDir = (lead(BmGeneStart) - BmGeneStart)) %>%
    summarise(sbum = sum(BmBlkDir, na.rm = T)) %>% 
    mutate(sbForward = ifelse(sbum <= 0, 0, 1)) %>% 
    select(sbid, sbForward)
  
  slimSyntRuns <- BmFuseSyntRuns %>%
    select(BmChr, FuseChr, sbid, FuseRunID, BmRunID, lastBmBlk, lastFuseBlk) %>% 
    group_by(sbid) %>% 
    slice(1) %>% 
    left_join(blockdir, by = 'sbid')
  
  ## Call blocks by their location in Bm and Fuse genomes
  slimSyntRuns <- slimSyntRuns %>%
    mutate(BmPos = ifelse(BmRunID == 1, 'b',
                          ifelse(BmRunID == lastBmBlk, 'e', 'i'))) %>% 
    mutate(FusePos = ifelse(FuseRunID == 1, 'b',
                          ifelse(FuseRunID == lastFuseBlk, 'e', 'i')))
  
  ## Encode end of chromosome matching
  slimSyntRuns <- slimSyntRuns %>% 
    mutate(matchEnd = ifelse(BmPos == 'i' | FusePos == 'i', 0,
                             ifelse(FusePos == 'b' & BmPos == 'b' & sbForward == 1, 1,
                                    ifelse(FusePos == 'e' & BmPos == 'e' & sbForward == 1, 1,
                                           ifelse(FusePos == 'b' & BmPos == 'e' & sbForward == 0, 1,
                                                  ifelse(FusePos == 'e' & BmPos == 'b' & sbForward == 0, 1,0))))))
}
endCounts <- vector()
for (i in 1:10000){
  fusedsynt.df <- reordChrom2(BmoriPnapiSynt)
  endpositions <- countEnds(fusedsynt.df)
  endCounts <- append(endCounts, sum(endpositions$matchEnd))
}

endCounts <- data_frame(endCounts)
pval <- endCounts %>% filter(endCounts >= 18) %>% summarise(pval = length(endCounts)/10000)
permEndPlot <- ggplot(endCounts) +
  geom_histogram(aes(x = endCounts), bins = 20) +
  geom_vline(aes(xintercept = ObservedEnds), color = 'red') +
  labs(x = "Correctly placed and oriented ends", y = '') +
#  theme_minimal() +
  theme(text = element_text(size=16),
        axis.title.x = element_text(size = 16)) +
  annotate(geom = 'text', x = 17, y = 1350, label = '18', color = 'red')

png(filename = 'syntBlockOrderC.png', width = 475, height = 335)
p
dev.off()
