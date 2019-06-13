## Random genome GO clustering
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
  return(data_frame(genomeNum = genomeNum, block = i, test1 = test1,
                    test2 = test2, test3 = test3, test4 = test4, size = size))
}

## Possibly parallel implementation of GO term enrichment comparison
library(parallel)
# Calculate the number of cores
no_cores <- detectCores()
# Initiate cluster
cl <- makeCluster(60)
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(topGO))
clusterExport(cl, c('expectedBlocks', 'PnapGenes', 'geneNames', 'geneID2GO'))
i <- 1
j <- 990
expectedGOs1 <- parLapply(cl, i:j, sigGOs)
proc.time()
i <- i+990
j <- j+990
expectedGOs2 <- parLapply(cl, i:j, sigGOs)
proc.time()
i <- i+990
j <- j+990
expectedGOs3 <- parLapply(cl, i:j, sigGOs)
proc.time()
i <- i+990
j <- j+990
expectedGOs4 <- parLapply(cl, i:j, sigGOs)
proc.time()
i <- i+990
j <- j+990
expectedGOs5 <- parLapply(cl, i:j, sigGOs)
proc.time()
i <- i+990
j <- j+990
expectedGOs6 <- parLapply(cl, i:j, sigGOs)
proc.time()
i <- i+990
j <- j+990
expectedGOs7 <- parLapply(cl, i:j, sigGOs)
proc.time()
i <- i+990
j <- j+990
expectedGOs8 <- parLapply(cl, i:j, sigGOs)
proc.time()
i <- i+990
j <- j+990
expectedGOs9 <- parLapply(cl, i:j, sigGOs)
proc.time()
i <- i+990
j <- j+990
expectedGOs10 <- parLapply(cl, i:j, sigGOs)
proc.time()
stopCluster(cl)
