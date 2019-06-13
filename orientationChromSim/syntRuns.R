## Remove synt genes that occur in <6 consecutive length runs

library(dplyr)
library(data.table)

# Read data from blast searches
BmoriPnapiSynt <- fread('BmoriPnapiSyntALL_chrOnly.csv',
                        col.names = c('BmGene','BmChr','BmGeneStart','BmGeneEnd',
                                      'PnGene','PnChr','PnGeneStart','PnGeneEnd'))

# Sort synt Blocks by Pnapi_chr and PnStart
# Itteratively count runs of syntenic genes, remove runs of 1:5, print number of blocks

for (rl in 1:5){
  runlengths <- BmoriPnapiSynt %>% arrange(PnChr, PnGeneStart) %>% 
    mutate(uBmChr = paste(BmChr, PnChr, sep=".")) %>%
    .$uBmChr %>% as.vector() %>% rle()
  BmoriPnapiSynt <- BmoriPnapiSynt %>% arrange(PnChr, PnGeneStart) %>%
    mutate(runlen = rep(runlengths$lengths, runlengths$lengths)) %>%
    mutate(runid = rleid(BmChr)) %>% filter(runlen > rl)
#  print(max(BmoriPnapiSynt$runid))
}

# Reassign new run lengths and IDs by PnChr
runlengths <- BmoriPnapiSynt %>% arrange(PnChr, PnGeneStart) %>% 
  .$BmChr %>% as.vector() %>% rle()
BmoriPnapiSynt <- BmoriPnapiSynt %>% arrange(PnChr, PnGeneStart) %>% 
  mutate(runlen = rep(runlengths$lengths, runlengths$lengths)) %>%
  group_by(PnChr) %>% mutate(runid = rleid(BmChr))
BmoriPnapiSynt <- BmoriPnapiSynt %>% ungroup(PnChr) %>% mutate(sbid = rleid(BmChr)) %>% 
  group_by(BmChr) %>% arrange(BmGeneStart) %>% mutate(BmRunID = rleid(sbid)) %>% ungroup()

# Generage syntenic block numbers from Bmori chromosomal perspective
for (rl in 1:5){
  runlengths <- BmoriPnapiSynt %>% arrange(BmChr, BmGeneStart) %>% 
    mutate(uPnChr = paste(PnChr, BmChr, sep=".")) %>%
    .$uPnChr %>% as.vector() %>% rle()
  BmoriPnapiSyntX <- BmoriPnapiSynt %>% arrange(BmChr, BmGeneStart) %>%
    mutate(Bmrunlen = rep(runlengths$lengths, runlengths$lengths)) %>%
    mutate(Bmrunid = rleid(PnChr)) %>% filter(runlen > rl)
  #  print(max(BmoriPnapiSynt$runid))
}

runlengths <- BmoriPnapiSyntX %>% arrange(BmChr, BmGeneStart) %>% 
  .$PnChr %>% as.vector() %>% rle()
BmoriPnapiSyntX <- BmoriPnapiSyntX %>% arrange(BmChr, BmGeneStart) %>% 
  mutate(Bmrunlen = rep(runlengths$lengths, runlengths$lengths)) %>%
  group_by(BmChr) %>% mutate(Bmrunid = rleid(PnChr))
BmoriPnapiSyntX <- BmoriPnapiSyntX %>% ungroup(BmChr) %>% mutate(Bmsbid = rleid(PnChr)) %>% 
  group_by(PnChr) %>% arrange(PnGeneStart) %>% mutate(PnRunID = rleid(Bmsbid)) %>% ungroup()


BmoriPnapiSyntX %>% 
  select(BmGene, BmChr, BmGeneStart, Bmsbid, PnGene, PnChr, PnGeneStart, sbid) %>%
  arrange(BmChr, BmGeneStart) %>% 
  fwrite(file = 'geneSbID_BmoriPnapi.tab', sep = '\t')
