### Results calculations
library(data.table)
library(dplyr)

#average orthologs in blocks
BmoriPnapiSynt %>% arrange(PnChr, PnGeneStart) %>% group_by(sbid) %>% slice(1) %>% ungroup() %>% .$runlen %>% mean()
BmoriPnapiSynt %>% arrange(PnChr, PnGeneStart) %>% group_by(PnChr) %>% 
  summarise(numblocks = max(runid)) %>% summarise(avgblocks = mean(numblocks), sdblocks = sd(numblocks))
#average genes in blocks on napi
GeneBlocks %>% arrange(PnChr, PnGeneStart) %>% group_by(sbid) %>% 
  summarise(blocksize = n()) %>% summarise(avgsize = mean(blocksize), stddev = sd(blocksize))
GeneBlocks %>% arrange(PnChr, PnGeneStart) %>% group_by(sbid) %>% 
  summarise(size = max(PnGeneEnd) - min(PnGeneStart)) %>% 
  summarise(avgsize = mean(size), stddev = sd(size))

