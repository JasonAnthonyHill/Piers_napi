### Simulate the merging of Bmori syntenic blocks
### and examine for maintainence of orientation

### Visualize by overlaying scaffold and syntenic block figures
### with markers signifying genes that come from varying distances of
### the ends of B. mori chromosomes

library(ggplot2)
library(dplyr)
library(data.table)

### Read in syntenic gene location and chromosome lengths and merge
BmoriPnapiSynt <- fread('BmoriPnapiSyntALL_chrOnly.csv',
                        col.names = c('BmGene','BmChr','BmGeneStart','BmGeneEnd',
                                      'PnGene','PnChr','PnGeneStart','PnGeneEnd'))
cLBmori <- fread('Bmori_chr_lengths.txt',
                     col.names = c('BmChr','BmChrLen'))

cLPnapi <- fread('Pnapi_chr_lengths.txt',
                 col.names = c('PnChr','PnChrLen'))

BPS <- BmoriPnapiSynt %>% left_join(cLBmori, by = 'BmChr') %>% left_join(cLPnapi, by = 'PnChr')

### Look for genes close to ends of BM
BPS %>% filter(BmGeneStart < 50000 | BmGeneEnd > (BmChrLen - 50000))

## Appearantly they do exist will plot on scaffold, syntenic block overlay
source('custom_arrow.R')

dSB <- fread('directionalSyntBlocks.csv')
scaffs <- fread('ScaffChrom_id.txt')

# custPal <- fread(file = "BmChrColors.txt", header = F)
# mycolors <- custPal$V2
# names(mycolors) <- custPal$V1

### Distance from end gene map
BPS <- BPS %>% mutate(dfeBm = ifelse(BmGeneStart < BmChrLen/2, BmGeneStart, BmChrLen - BmGeneEnd)) %>% 
  mutate(dfePn = ifelse(PnGeneStart < PnChrLen/2, PnGeneStart, PnChrLen - PnGeneEnd)) %>% 
  mutate(pdfeBm = dfeBm/BmChrLen) %>% mutate(pdfePn = dfePn/PnChrLen)

ggplot(BPS) +
  geom_density2d(aes(dfePn, dfeBm), color = "black") +
  viridis::scale_fill_viridis() +
  theme_minimal()

ggplot(BPS) +
  stat_density_2d(aes(dfePn, dfeBm, fill = ..level..), geom = 'polygon') +
  viridis::scale_fill_viridis() +
  theme_minimal()

ggplot(BPS) +
  geom_bin2d(aes(dfePn, dfeBm), color = "black", bins = 10) +
  viridis::scale_fill_viridis() +
  theme_minimal() +
  labs(x = "% Distance from end - Pnapi", y = '% Distance from end - Bmori') +
  guides(fill=guide_legend(title="Genes"))

