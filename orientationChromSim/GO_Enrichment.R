## Extract gene ranges and GO terms for each syntenic block and check
## for enrichment

library(ggplot2)
library(plyr)
library(data.table)
library(tidyr)
library(dplyr)

source('syntRuns.R')
PnapGeneLoc <- fread('PnapGeneLocation.txt', 
                     col.names = c('PnGene', 'PnChr', 'PnGeneStart', 'PnGeneEnd'),
                     colClasses = c('character', 'character', 'integer', 'integer')) %>% 
  arrange(PnChr, PnGeneStart)
PnapGeneLoc$PnChr <- gsub('Chromosome', 'Pnapi_chr', PnapGeneLoc$PnChr)

PnapGOtoGene <- fread('rc3GO.tsv', col.names = c('PnGene', 'GO_term'))

blockStartEnd <- BmoriPnapiSynt %>% 
  arrange(sbid, PnGeneStart) %>%
  group_by(sbid) %>% 
  summarise(blockStart = first(PnGeneStart), blockEnd = last(PnGeneStart), PnChr = first(PnChr))

## Go through each sblock and extract the genes between the start and end

getBlocks.func <- function(x){
  sbid.f <- as.numeric(x[1])
  blockStart.f <- as.numeric(x[2])
  blockEnd.f <- as.numeric(x[3])
  PnChr.f <- x[4]
  block.df <- PnapGeneLoc %>%
    arrange(PnChr, PnGeneStart) %>% 
    filter(PnChr == PnChr.f) %>% 
    filter(PnGeneStart >= blockStart.f & PnGeneEnd <= blockEnd.f) %>% 
    mutate(sbid = sbid.f)
  return(block.df)
}

GeneBlocks <- apply(blockStartEnd, 1, getBlocks.func) %>% bind_rows()
GeneBlocks <- GeneBlocks %>% left_join(PnapGOtoGene, by = 'PnGene')

PnapGOtoGene %>% 
  group_by(PnGene) %>% 
  mutate(id = row_number()) %>% 
  spread(id, GO_term) %>%
  unite(val, 2:103, sep = ',') %>% 
  mutate(val = gsub(',NA', '', val)) %>% View()
  fwrite(file = 'gene2GO.txt', sep = '\t', col.names = F)
  
library(topGO)
geneID2GO <- readMappings(file = 'gene2GO.txt')
geneNames <- names(geneID2GO)
resultsTable <- data_frame()
for(i in 1:99){
myInterestingGenes <- GeneBlocks %>% 
  filter(sbid == i) %>% 
  dplyr::select(PnGene) %>%
  unlist()
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
GOdata <- new("topGOdata",
              description = base::paste("GO Topology of Syntenic Block",i, sep = ' '), 
              ontology = "BP",
              allGenes = geneList,
#              geneSel = myInterestingGenes,
              nodeSize = 5,
              annot = annFUN.gene2GO, 
              gene2GO = geneID2GO)

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
tableRes <- GenTable(GOdata, classicFisher = resultFisher, topNodes = 50)
tableRes <- tableRes %>% mutate(sbid = i)
resultsTable <- resultsTable %>% bind_rows(tableRes)
}

obsTest1 <- resultsTable %>% 
  filter(Significant >= 3 & classicFisher < 0.05) %>% 
  summarise(ndist = n_distinct(sbid))
obsTest2 <- resultsTable %>% 
  filter(Significant >= 3 & classicFisher < 0.01) %>% 
  summarise(ndist = n_distinct(sbid))
obsTest3 <- resultsTable %>% 
  filter(Significant >= 5 & classicFisher < 0.05) %>% 
  summarise(ndist = n_distinct(sbid))
obsTest4 <- resultsTable %>% 
  filter(Significant >= 5 & classicFisher < 0.01) %>% 
  summarise(ndist = n_distinct(sbid))
observedSigBlocks <- data_frame(test = c('test1bool','test2bool','test3bool','test4bool'), 
                                observedSignif = as.numeric(c(obsTest1, obsTest2, obsTest3, obsTest4)))
## Random genome GO clustering
# Assign genes/GO terms to each fraction
PnapGenes <- PnapGeneLoc %>% mutate(geneNumber = row_number())

# Generate large matrix of 10000 genomes composed of 99 random normally distributed fragments adding up to 1
jobFun <- function(n) {
  m <- matrix(runif(99*n,0.25,.75), ncol=99)
  m<- sweep(m, 1, rowSums(m), FUN="/")
  m
}
JOB <- jobFun(10000)
eblocks <- JOB*13622
# Generate a matrix of genomes composed of 99 blocks drawn from observed distribution
# pull <- function(x,y) {x[,if(is.name(substitute(y))) deparse(substitute(y)) else y, drop = FALSE][[1]]}
bsvec <- GeneBlocks %>% arrange(sbid) %>% .$sbid %>% rle() %>% .$lengths
job2Fun <- function(x){
  m <- replicate(x, base::sample(bsvec,99))
  m %>% t() %>% return()
}
eblocks <- job2Fun(10000)

## Calculate GO enrichment in expected blocks
expectedBlocks <- eblocks %>% as.data.frame()
colnames(expectedBlocks) <- paste(seq(1:99))
expectedBlocks <- expectedBlocks %>%
  mutate(genomeNum = row_number()) %>% 
  gather(sbid, size, 1:99) %>%
  mutate(sbid = as.numeric(sbid)) %>% 
  mutate(size = round(size)) %>% 
  arrange(genomeNum, sbid) %>% 
  group_by(genomeNum) %>% 
  mutate(lastGene = cumsum(size)) %>% 
  mutate(firstGene = (lastGene - size) + 1)

sigGOs <- function(i){
  first <- expectedBlocks$firstGene[i]
  last <- expectedBlocks$lastGene[i]
  genomeNum <- expectedBlocks$genomeNum[i]
  # print(first)
myInterestingGenes <- PnapGenes %>% 
  filter(geneNumber >= first & geneNumber <= last) %>% 
  dplyr::select(PnGene) %>%
  unlist()
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
base::levels(geneList) <- c('0','1')
GOdata <- new("topGOdata",
              description = base::paste("GO Topology of Syntenic Block",i, sep = ' '), 
              ontology = "BP",
              allGenes = geneList,
              #              geneSel = myInterestingGenes,
              nodeSize = 5,
              annot = annFUN.gene2GO, 
              gene2GO = geneID2GO)

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
tableRes <- GenTable(GOdata, classicFisher = resultFisher, topNodes = 50)
# tableRes <- tableRes %>% mutate(sbid = i)
# resultsTable <- resultsTable %>% bind_rows(tableRes)

test1 <- tableRes %>% 
  filter(Significant >= 3 & classicFisher < 0.05) %>% 
  summarise(ndist = n_distinct(GO.ID)) %>% 
  .$ndist

test2 <- tableRes %>% 
  filter(Significant >= 3 & classicFisher < 0.01) %>% 
  summarise(ndist = n_distinct(GO.ID)) %>% 
  .$ndist

test3 <- tableRes %>% 
  filter(Significant >= 5 & classicFisher < 0.05) %>% 
  summarise(ndist = n_distinct(GO.ID)) %>% 
  .$ndist

test4 <- tableRes %>% 
  filter(Significant >= 5 & classicFisher < 0.01) %>% 
  summarise(ndist = n_distinct(GO.ID)) %>% 
  .$ndist

size <- last - first
print(i)
return(data_frame(genomeNum = genomeNum, block = i, test1 = test1, 
                  test2 = test2, test3 = test3, test4 = test4, size = size))
}

## Possibly parallel implementation of GO term enrichment comparison
library(parallel)
# Calculate the number of cores
no_cores <- detectCores()
# Initiate cluster
cl <- makeCluster(7)
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(topGO))
clusterExport(cl, c('expectedBlocks', 'PnapGenes', 'geneNames', 'geneID2GO'))

timestart <- proc.time()
expectedGOs2 <- parLapply(cl, 1:nrow(expectedBlocks), sigGOs)
timestop <- proc.time()

stopCluster(cl)

# Single cor working version

expectedGOs <- expectedBlocks %>%
  ungroup() %>%
  group_by(genomeNum, sbid) %>%
  do(sigGOs(.$firstGene,.$lastGene))
#bind_rows(expectedGOs, expectedGOs2) %>% fwrite(file = 'expectedGOs.csv')

# load('expected5000.R')
load('expected10k_data.R')
expectedGOs <- expectedGOs10k %>% bind_rows()


eGOs <- expectedGOs %>% 
  mutate(test1bool = ifelse(test1 > 0, 1, 0)) %>%
  mutate(test2bool = ifelse(test2 > 0, 1, 0)) %>% 
  mutate(test3bool = ifelse(test3 > 0, 1, 0)) %>% 
  mutate(test4bool = ifelse(test4 > 0, 1, 0))

sumGOs <- eGOs %>% 
  group_by(genomeNum) %>% 
  summarise_each(funs(sum),test1bool,test2bool,test3bool,test4bool) %>% 
  gather('test','expectedSignif',2:5) %>% 
  mutate(observedSignif = ifelse(test == 'test1bool', obsTest1,
                                 ifelse(test == 'test2bool', obsTest2,
                                        ifelse (test == 'test3bool', obsTest3, obsTest4))))

ggplot(sumGOs) + 
  geom_histogram(aes(x = expectedSignif), bins = 20, alpha = .75) +
  geom_vline(aes(xintercept = observedSignif), data = observedSigBlocks, color = 'red') +
  coord_cartesian() +
  theme_minimal() +
  labs(x = 'Signif GO enriched blocks') +
  facet_wrap(~test)

txtsize <- 84
permGOPlot <- sumGOs %>% 
  filter(test == 'test2bool') %>% 
ggplot() + 
  geom_histogram(aes(x = expectedSignif), 
                 binwidth = 1,
                 fill = "black",
                 col = "white",
                 size = 2) +
  geom_segment(x = 57, y = 0, xend= 57, yend = 600, color = 'red', size = 3) +
  coord_cartesian(xlim = c(25,60), expand = T) +
  # theme_minimal() +
  theme_bw() +
  theme(text = element_text(size=txtsize, face = "bold"),
        axis.title.x = element_text(size = txtsize, face="bold", margin = margin(t = txtsize*.5), color = "gray25"),
        axis.title.y = element_text(size = txtsize, face="bold", margin = margin(r = txtsize*.5), color = "gray25"),
        axis.line = element_line(colour = 'black', size = 3)) +
  labs(x = 'GO term enriched blocks', y = 'Occurences in 10k Sims.') +
  annotate(geom = 'text', x = 55, y = 500, label = '57', color = 'red', size = txtsize*.30) +
  expand_limits(x = 0, y = 0)
png(filename = 'GOfig.png', width = 1600, height = 1280)
permGOPlot
dev.off()


#Results table
GOresultsTable <- resultsTable %>% 
  filter(Significant >= 3, classicFisher < 0.01) %>%
  left_join(blockStartEnd)
write.table(GOresultsTable, 'GOresultsTable.tsv', sep = '\t', row.names = F)
