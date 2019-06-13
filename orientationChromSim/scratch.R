ggplot(condensed.df) + 
  # geom_histogram(aes(fuseGene.dfe/fuseGeneNumMax))
  # geom_histogram(aes(PnGene.dfe/PnGeneNumMax))
  geom_histogram(aes(BmGene.dfe/BmGeneNumMax))

y.df <- condensed.df %>%
  ggplot() +
  geom_bin2d(aes(BmGene.dfe/BmGeneNumMax, fuseGene.dfe/fuseGeneNumMax), color = "black", bins = 25, na.rm = T) +
  # geom_bin2d(aes(BmGene.dfe, fuseGene.dfe), color = "black", bins = 25) +
  viridis::scale_fill_viridis() +
  theme_minimal() +
  labs(x = "% Distance from end - fuse", y = '% Distance from end - Bmori') +
  guides(fill=guide_legend(title="Genes")) +
  coord_fixed()

x.df <- ggplot_build(y.df)

base::lapply(expectedBlocks, sigGOs(expectedBlocks$firstGene, expectedBlocks$lastGene))
base::lapply(expectedBlocks, print(expectedBlocks$firstGene))

no_cores <- detectCores() - 2
cl <- makeCluster(no_cores)
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(topGO))
clusterExport(cl, c('expectedBlocks', 'PnapGenes', 'geneNames', 'geneID2GO'))

expectedGOs2 <- parLapply(cl, 1:nrow(expectedBlocks), sigGOs(i))
stopCluster(cl)

no_cores <- detectCores() - 2
cl <- makeCluster(no_cores)
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(topGO))
testdf <- data_frame(a = c(1,2,3), b = c(4,5,6), c = c(7,8,9))
clusterExport(cl, 'testdf')
testfun <- function(i){
  a <- testdf$a[i]
  b <- testdf$b[i]
  ab <- data_frame(a = a, b = b)
  return(ab)
}
testmulti <- parLapply(cl, 1:nrow(testdf), testfun)
lapply(1:nrow(testdf), testfun)
stopCluster(cl)

BmoriPnapiSynt %>% arrange(PnChr, PnGeneStart) %>%  group_by(PnChr, BmChr) %>% filter(row_number() == 1) %>% View()

## Load all synt and look for other genes on Bombyx Y
AllOrthoSyntPnBm <- fread('BmoriPnapiSyntALL.csv', col.names = c('BmGene','BmChr','BmStart','BmEnd',
                                                                 'PnGene','PnChr','PnStart','PnEnd'))

AllOrthoSyntPnBm %>% arrange(BmChr,BmStart) %>% filter(BmChr == 'Bmori_chr_1') %>% count(PnChr)

MPscaffJoins <- fread('../ReadSpanSynt/MPspanScaffJoins.tsv')
MPscaffJoins %>% filter(V7 > 0) %>% summarise(avg = mean(V7), vari = sqrt(var(V7)))
ggplot(MPscaffJoins) + geom_line(aes(x = row_number(V7), y = V7))

BmoriPnapiSyntx <- fread('BmoriPnapiSyntALL_chrOnly.csv',
                        col.names = c('BmGene','BmChr','BmGeneStart','BmGeneEnd',
                                      'PnGene','PnChr','PnGeneStart','PnGeneEnd'))
# Reassign new run lengths and IDs by PnChr
runlengthsx <- BmoriPnapiSyntx %>% arrange(PnChr, PnGeneStart) %>% 
  .$BmChr %>% as.vector() %>% rle()
BmoriPnapiSyntx <- BmoriPnapiSyntx %>% arrange(PnChr, PnGeneStart) %>% 
  mutate(runlen = rep(runlengthsx$lengths, runlengthsx$lengths)) %>%
  group_by(PnChr) %>% mutate(runid = rleid(BmChr))
BmoriPnapiSyntx <- BmoriPnapiSyntx %>% ungroup(PnChr) %>% mutate(sbid = rleid(BmChr)) %>% 
  group_by(BmChr) %>% arrange(BmGeneStart) %>% mutate(BmRunID = rleid(sbid)) %>% ungroup()
BmoriPnapiSyntx %>% summarise(total = max(sbid))
BmoriPnapiSynt %>% summarise(total = max(sbid))
