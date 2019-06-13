## Generate null gene to end distribution with random syntenic joins

library(dplyr)
library(data.table)
library(ggplot2)

source('syntRuns.R')

## Assign gene order numbers to Bm and Pn genomes

BmoriPnapiSynt <- BmoriPnapiSynt %>%
  arrange(BmChr, BmGeneStart) %>%
  group_by(BmChr) %>% 
  mutate(BmGeneNum = row_number()) %>% 
  ungroup() %>% 
  arrange(PnChr, PnGeneStart) %>% 
  group_by(PnChr) %>% 
  mutate(PnGeneNum = row_number())

## Function to rearrange syneny data frame by random draw
reordChrom <- function(df.in){
#  df.in <- BmoriPnapiSynt
  sbids <- df.in$sbid %>% unique()
  ord <- sample(sbids, replace = FALSE, length(sbids))
  gn.df <- split(ord, ceiling(seq_along(ord)/4)) %>% 
    plyr::ldply(data.frame, .id = 'fuseChr') %>% 
    rename(sbid = X..i..) %>% 
    mutate(orient = sample(c('f','r'), length(sbids), replace = T)) %>% 
    right_join(df.in, by = 'sbid') %>% 
    mutate(sbid = factor(sbid, levels = ord)) %>% 
    group_by(sbid) %>%  
    arrange(sbid, ifelse(orient == "f", BmGeneStart, -BmGeneStart)) %>%
    ungroup(sbid) %>%
    group_by(fuseChr) %>% 
    mutate(fuseGeneNum = row_number()) %>%
    ungroup()
  
  maxGnBm.df <- gn.df %>% 
    group_by(BmChr) %>% 
    slice(which.max(BmGeneNum)) %>% 
    select(BmChr, BmGeneNum) %>% 
    arrange(BmChr) %>% 
    ungroup() %>% 
    rename(BmGeneNumMax = BmGeneNum)
  
  maxGnPn.df <- gn.df %>% 
    group_by(PnChr) %>% 
    slice(which.max(PnGeneNum)) %>% 
    select(PnChr, PnGeneNum) %>% 
    arrange(PnChr) %>% 
    ungroup() %>% 
    rename(PnGeneNumMax = PnGeneNum)
  
  maxGnFuse.df <- gn.df %>% 
    group_by(fuseChr) %>% 
    slice(which.max(fuseGeneNum)) %>% 
    select(fuseChr, fuseGeneNum) %>% 
    arrange(fuseChr) %>% 
    ungroup() %>% 
    rename(fuseGeneNumMax = fuseGeneNum)
  
  condensed.df <- gn.df %>% left_join(maxGnBm.df, by = 'BmChr') %>% 
    left_join(maxGnPn.df, by = 'PnChr') %>% 
    left_join(maxGnFuse.df, by = 'fuseChr') %>% 
    mutate(BmGene.dfe = ifelse(BmGeneNum < BmGeneNumMax/2,
                               BmGeneNum, (BmGeneNumMax - BmGeneNum))) %>% 
    mutate(PnGene.dfe = ifelse(PnGeneNum < PnGeneNumMax/2,
                               PnGeneNum, (PnGeneNumMax - PnGeneNum))) %>% 
    mutate(fuseGene.dfe = ifelse(fuseGeneNum < fuseGeneNumMax/2,
                                 fuseGeneNum, (fuseGeneNumMax - fuseGeneNum)))# %>% 
    # select(BmChr, BmGeneNum, BmGeneNumMax, BmGene.dfe,
    #        fuseChr, fuseGeneNum, fuseGeneNumMax, fuseGene.dfe)
  
  return(condensed.df)
}

## Make list of 1000 Estimated distance tables
## percent by order 25 bins
Epo25 <- data_frame()
for (i in 1:10000){
  test <- reordChrom(BmoriPnapiSynt) %>% arrange(fuseChr, fuseGeneNum)
  p <-test %>% filter() %>% ggplot() +
    geom_bin2d(aes(BmGene.dfe/BmGeneNumMax, fuseGene.dfe/fuseGeneNumMax), color = "black", bins = 25) +
    viridis::scale_fill_viridis() +
    theme_minimal() +
    labs(x = "% Distance from end - fuse", y = '% Distance from end - Bmori') +
    guides(fill=guide_legend(title="Genes")) +
    coord_fixed()
  q <- ggplot_build(p)
  Epo25 <- q$data %>% data.frame() %>% select(xbin, ybin, count) %>% 
    rename(Estimated = count, x = xbin, y = ybin) %>% bind_rows(Epo25)
}

## percent by order 10 bins
Epo10 <- data_frame()
for (i in 1:100){
  test <- reordChrom(BmoriPnapiSynt) %>%
    arrange(fuseChr, fuseGeneNum)
  p <-test %>% filter() %>% ggplot() +
    geom_bin2d(aes(BmGene.dfe/BmGeneNumMax, fuseGene.dfe/fuseGeneNumMax), color = "black", bins = 10) +
    viridis::scale_fill_viridis() +
    theme_minimal() +
    labs(x = "% Distance from end - fuse", y = '% Distance from end - Bmori') +
    guides(fill=guide_legend(title="Genes")) +
    coord_fixed()
  q <- ggplot_build(p)
  Epo10 <- q$data %>% data.frame() %>% select(xbin, ybin, count) %>% 
    rename(Estimated = count, x = xbin, y = ybin) %>% bind_rows(Epo10)
}

## order 10 bins
Eo10 <- data_frame()
for (i in 1:10000){
  test <- reordChrom(BmoriPnapiSynt) %>% arrange(fuseChr, fuseGeneNum)
  p <-test %>% filter() %>% ggplot() +
    geom_bin2d(aes(BmGene.dfe, fuseGene.dfe), color = "black", bins = 10) +
    viridis::scale_fill_viridis() +
    theme_minimal() +
    labs(x = "% Distance from end - fuse", y = '% Distance from end - Bmori') +
    guides(fill=guide_legend(title="Genes")) +
    coord_fixed()
  q <- ggplot_build(p)
  Eo10 <- q$data %>% data.frame() %>% select(xbin, ybin, count) %>% 
    rename(Estimated = count, x = xbin, y = ybin) %>% bind_rows(Eo10)
}

## order 25 bins
Eo25 <- data_frame()
for (i in 1:10000){
  test <- reordChrom(BmoriPnapiSynt) %>% arrange(fuseChr, fuseGeneNum)
  p <-test %>% filter() %>% ggplot() +
    geom_bin2d(aes(BmGene.dfe, fuseGene.dfe), color = "black", bins = 25) +
    viridis::scale_fill_viridis() +
    theme_minimal() +
    labs(x = "% Distance from end - fuse", y = '% Distance from end - Bmori') +
    guides(fill=guide_legend(title="Genes")) +
    coord_fixed()
  q <- ggplot_build(p)
  Eo25 <- q$data %>% data.frame() %>% select(xbin, ybin, count) %>% 
    rename(Estimated = count, x = xbin, y = ybin) %>% bind_rows(Eo25)
}

## Generate Observed data set for different bins and normalization
gn.df <- BmoriPnapiSynt %>% group_by(PnChr) %>% arrange(PnChr, PnGeneStart) %>% 
  mutate(PnGeneNum = row_number()) %>%
  ungroup(PnChr) %>% group_by(BmChr) %>% arrange(BmChr, BmGeneStart) %>% 
  mutate(BmGeneNum = row_number()) %>% ungroup()

maxGnBm.df <- gn.df %>% group_by(BmChr) %>% slice(which.max(BmGeneNum)) %>% 
  select(BmChr, BmGeneNum) %>% arrange(BmChr) %>% ungroup() %>% rename(BmGeneNumMax = BmGeneNum)

maxGnPn.df <- gn.df %>% group_by(PnChr) %>% slice(which.max(PnGeneNum)) %>% 
  select(PnChr, PnGeneNum) %>% arrange(PnChr) %>% ungroup() %>% rename(PnGeneNumMax = PnGeneNum)

condensed.df <- gn.df %>% left_join(maxGnBm.df, by = 'BmChr') %>% 
  left_join(maxGnPn.df, by = 'PnChr') %>% 
  mutate(BmGene.dfe = ifelse(BmGeneNum < BmGeneNumMax/2,
                             BmGeneNum, (BmGeneNumMax - BmGeneNum))) %>% 
  mutate(PnGene.dfe = ifelse(PnGeneNum < PnGeneNumMax/2,
                             PnGeneNum, (PnGeneNumMax - PnGeneNum))) %>% 
  select(BmChr, BmGeneStart, BmGeneNum, BmGeneNumMax, BmGene.dfe,
         PnChr, PnGeneStart, PnGeneNum, PnGeneNumMax, PnGene.dfe)

## Observed percent order 25 bins
Opo25.plot <- ggplot(condensed.df) +
  geom_bin2d(aes(BmGene.dfe/BmGeneNumMax, PnGene.dfe/PnGeneNumMax), color = "black", bins = 25) +
  viridis::scale_fill_viridis() +
  theme_minimal() +
  labs(x = "% Distance from end - Pnapi", y = '% Distance from end - Bmori') +
  guides(fill=guide_legend(title="Genes"))

Opo25.plot.build <- ggplot_build(Opo25.plot)

Opo25 <- Opo25.plot.build$data %>% data.frame() %>% select(xbin, ybin, count) %>% 
  rename(Observed = count, x = xbin, y = ybin)

## Observed percent order 10 bins
Opo10.plot <- ggplot(condensed.df) +
  geom_bin2d(aes(BmGene.dfe/BmGeneNumMax, PnGene.dfe/PnGeneNumMax), 
             color = "black", bins = 10) +
  viridis::scale_fill_viridis(limit = c(0, 200),
                              breaks = c(0,50,100,150,200)) +
  coord_fixed() +
  labs(x = "% DFE - Pnapi", y = '% DFE - Bmori') +
  guides(fill=guide_legend(title="Genes")) +
  scale_x_continuous(breaks = c(0,.1,.2,.3,.4,.5),
                     labels = c(0,10,20,30,40,50),
                     expand = c(0,0)) +
  scale_y_continuous(breaks = c(0,.1,.2,.3,.4,.5),
                     labels = c(0,10,20,30,40,50),
                     expand = c(0,0)) +
  labs(x = expression(paste('% DFE - ', italic('P.napi'))), 
       y = expression(paste('% DFE - ', italic('B.mori')))) +
  coord_fixed() +
  theme(text = element_text(size=16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text=element_text(size=12),
        legend.title = element_text(size=14))

Opo10.plot.build <- ggplot_build(Opo10.plot)

Opo10 <- Opo10.plot.build$data %>% data.frame() %>% select(xbin, ybin, count) %>% 
  rename(Observed = count, x = xbin, y = ybin)

## Observed order 25 bins
Oo25.plot <- ggplot(condensed.df) +
  geom_bin2d(aes(BmGene.dfe, PnGene.dfe), color = "black", bins = 25) +
  viridis::scale_fill_viridis() +
  theme_minimal() +
  labs(x = "% Distance from end - Pnapi", y = '% Distance from end - Bmori') +
  guides(fill=guide_legend(title="Genes"))

Oo25.plot.build <- ggplot_build(Oo25.plot)

Oo25 <- Oo25.plot.build$data %>% data.frame() %>% select(xbin, ybin, count) %>% 
  rename(Observed = count, x = xbin, y = ybin)

## Observed order 10 bins
Oo10.plot <- ggplot(condensed.df) +
  geom_bin2d(aes(BmGene.dfe, PnGene.dfe), color = "black", bins = 10) +
  viridis::scale_fill_viridis() +
  theme_minimal() +
  labs(x = "% Distance from end - Pnapi", y = '% Distance from end - Bmori') +
  guides(fill=guide_legend(title="Genes"))

Oo10.plot.build <- ggplot_build(Oo10.plot)

Oo10 <- Oo10.plot.build$data %>% data.frame() %>% select(xbin, ybin, count) %>% 
  rename(Observed = count, x = xbin, y = ybin)

## Combine estimated and observed data frames
EOpo25 <- left_join(Epo25, Opo25, by = c('x','y')) %>% group_by(x,y)
EOpo10 <- left_join(Epo10, Opo10, by = c('x','y')) %>% group_by(x,y)
EOo25 <- left_join(Eo25, Oo25, by = c('x','y')) %>% group_by(x,y)
EOo10 <- left_join(Eo10, Oo10, by = c('x','y')) %>% group_by(x,y)

## Save in list and write
# EOlist <- list(EOpo10, EOpo25, EOo10, EOo25)
# save(EOlist, file = "EOlist.data.R")
# Elist <- list(Epo10, Epo25, Eo10, Eo25)
# save(Elist, file = 'Elist.data.R')
## Load list of data frames and unpack to globalenv
EOlist <- load(file = "EOlist.data.R")
list2env(EOlist, .GlobalEnv)

## Derive P values from E-O comparisons result
resultEOo10 <- summarize(EOo10, p_left = sum(Observed >= Estimated)/n(),
                    p_right = sum(Observed <= Estimated)/n(),
                    Emean = mean(Estimated)) %>%
  mutate(p = pmin(p_left,p_right)*2,
         direction = p_left < p_right,
         dir_p = if_else(direction, p, -p))
resultEOo10$dir_p2 <- ifelse(resultEOo10$direction, 1 - resultEOo10$dir_p, -1 - resultEOo10$dir_p)

resultEOo25 <- summarize(EOo25, p_left = sum(Observed >= Estimated)/n(),
                         p_right = sum(Observed <= Estimated)/n(),
                         Emean = mean(Estimated)) %>%
  mutate(p = pmin(p_left,p_right)*2,
         direction = p_left < p_right,
         dir_p = if_else(direction, p, -p))
resultEOo25$dir_p2 <- ifelse(resultEOo25$direction, 1 - resultEOo25$dir_p, -1 - resultEOo25$dir_p)

resultEOpo10 <- summarize(EOpo10, p_left = sum(Observed >= Estimated)/n(),
                         p_right = sum(Observed <= Estimated)/n(),
                         Emean = mean(Estimated)) %>%
  mutate(p = pmin(p_left,p_right)*2,
         direction = p_left < p_right,
         dir_p = if_else(direction, p, -p))
resultEOpo10$dir_p2 <- ifelse(resultEOpo10$direction, 1 - resultEOpo10$dir_p, -1 - resultEOpo10$dir_p)

resultEOpo25 <- summarize(EOpo25, p_left = sum(Observed >= Estimated)/n(),
                         p_right = sum(Observed <= Estimated)/n(),
                         Emean = mean(Estimated)) %>%
  mutate(p = pmin(p_left,p_right)*2,
         direction = p_left < p_right,
         dir_p = if_else(direction, p, -p))
resultEOpo25$dir_p2 <- ifelse(resultEOpo25$direction, 1 - resultEOpo25$dir_p, -1 - resultEOpo25$dir_p)

## Plot results

cols <- RColorBrewer::brewer.pal(10, "PuOr")
vals <- c(1, .95, .9, .8, 0)
vals <- c(-vals, rev(vals))
permEopoPlot <- ggplot(resultEOpo10) +
  geom_tile(aes(x, y, fill = dir_p2)) +
  geom_point(data = . %>% filter(dir_p2 <= -0.95 | dir_p2 >= 0.95),
             aes(x,y), color = 'white', fill = 'white', shape = 18, show.legend = FALSE) +
  scale_fill_gradientn(colours = cols, values = scales::rescale(vals), limit = c(-1, 1),
                       guide = guide_colorbar(barheight = 15),
                       breaks = c(-1, -.95, -.9, -.8, 0, .8, .9, .95, 1),
                       labels = c(0, 0.05, 0.1, 0.2, 0, 0.2, 0.1, 0.05, 0)) +
  scale_x_continuous(breaks = c(0,2,4,6,8,10),
                     labels = c(0,10,20,30,40,50),
                     expand = c(0,0)) +
  scale_y_continuous(breaks = c(0,2,4,6,8,10),
                     labels = c(0,10,20,30,40,50),
                     expand = c(0,0)) +
  labs(x = expression(paste('% DFE - ', italic('P.napi'))), 
       y = expression(paste('% DFE - ', italic('B.mori')))) +
  guides(fill=guide_legend(title = 'p-val')) +
  coord_fixed() +
  theme(text = element_text(size=16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text=element_text(size=12),
        legend.title = element_text(size=14))

## Generate Observed data set
gn.df <- BmoriPnapiSynt %>% group_by(PnChr) %>% arrange(PnChr, PnGeneStart) %>% 
  mutate(PnGeneNum = row_number()) %>%
  ungroup(PnChr) %>% group_by(BmChr) %>% arrange(BmChr, BmGeneStart) %>% 
  mutate(BmGeneNum = row_number()) %>% ungroup()

maxGnBm.df <- gn.df %>% group_by(BmChr) %>% slice(which.max(BmGeneNum)) %>% 
  select(BmChr, BmGeneNum) %>% arrange(BmChr) %>% ungroup() %>% rename(BmGeneNumMax = BmGeneNum)

maxGnPn.df <- gn.df %>% group_by(PnChr) %>% slice(which.max(PnGeneNum)) %>% 
  select(PnChr, PnGeneNum) %>% arrange(PnChr) %>% ungroup() %>% rename(PnGeneNumMax = PnGeneNum)

condensed.df <- gn.df %>% left_join(maxGnBm.df, by = 'BmChr') %>% 
  left_join(maxGnPn.df, by = 'PnChr') %>% 
  mutate(BmGene.dfe = ifelse(BmGeneNum < BmGeneNumMax/2,
                             BmGeneNum, (BmGeneNumMax - BmGeneNum))) %>% 
  mutate(PnGene.dfe = ifelse(PnGeneNum < PnGeneNumMax/2,
                             PnGeneNum, (PnGeneNumMax - PnGeneNum))) %>% 
  select(BmChr, BmGeneStart, BmGeneNum, BmGeneNumMax, BmGene.dfe,
         PnChr, PnGeneStart, PnGeneNum, PnGeneNumMax, PnGene.dfe)

p1 <- ggplot(condensed.df) +
  geom_bin2d(aes(BmGene.dfe, PnGene.dfe), color = "black", bins = 1) +
  viridis::scale_fill_viridis() +
  theme_minimal() +
  labs(x = "Order distance from end - Pnapi", y = 'Order distance from end - Bmori') +
  guides(fill=guide_legend(title="Genes"))

q1 <- ggplot_build(p1)

O <- q1$data %>% data.frame() %>% select(xbin, ybin, count) %>% 
  rename(Observed = count, x = xbin, y = ybin)

## Combine estimated and observed plots
EO <- left_join(E, O, by = c('x','y')) %>% group_by(x,y)

## Summarise result
result <- summarize(EO, p_left = sum(Observed <= Estimated)/n(),
                    p_right = sum(Observed >= Estimated)/n(),
                    Emean = mean(Estimated)) %>%
  mutate(p = pmin(p_left,p_right)*2,
         direction = p_left < p_right,
         dir_p = if_else(direction, p, -p))

cols <- RColorBrewer::brewer.pal(10, "PuOr")
vals <- c(1, .95, .9, .8, .5)
vals <- c(-vals, rev(vals))
result$dir_p2 <- ifelse(result$direction, 1 - result$dir_p, -1 - result$dir_p)
ggplot(result) +
  geom_tile(aes(x, y, fill = dir_p2)) +
  scale_fill_gradientn(colours = cols, values = scales::rescale(vals), limit = c(-1, 1),
                       guide = guide_colorbar(barheight = 15),
                       breaks = c(-1, -.9, -.5, 0, .5, .9, 1),
                       labels = c(0, 0.1, 0.5, 1, 0.5, 0.1, 0))



### Analysis using genomic coordinates to asess distance from the end of the chromosome



ggplot(result) +
  geom_tile(aes(x,y, fill = Emean)) +
  viridis::scale_fill_viridis()

test %>% filter(fuseChr ==1) %>% ggplot() +
  geom_point(aes(BmGene.dfe/BmGeneNumMax, fuseGene.dfe/fuseGeneNumMax, color = BmChr),
             size = .5)

test %>% filter(BmGene.dfe/BmGeneNumMax < 0.01 & fuseGene.dfe/fuseGeneNumMax < 0.01)

ggplot(test) +
  geom_histogram(aes(x = fuseGene.dfe/fuseGeneNumMax), bins = 50)
ggplot(test) +
  geom_histogram(aes(x = BmGene.dfe/BmGeneNumMax), bins = 50)

test <- reordChrom(BmoriPnapiSynt) %>% arrange(fuseChr, fuseGeneNum)
ggplot(test) +
  geom_point(aes(BmGene.dfe/BmGeneNumMax, fuseGene.dfe/fuseGeneNumMax),
             size = .5)

test <- chr %>% filter(newChrNum == 1) %>%
  group_by(sbid) %>%  
  arrange(sbid, ifelse(orient == "f", BmGeneStart, -BmGeneStart)) 

blo <- BmoriPnapiSynt %>%
  mutate(rosb =  factor(sbid, levels = ord)) %>%
  arrange(rosb)

gn.df <- BmoriPnapiSynt %>% group_by(PnChr) %>% arrange(PnChr, PnGeneStart) %>% 
  mutate(PnGeneNum = row_number()) %>%
  ungroup(PnChr) %>% group_by(BmChr) %>% arrange(BmChr, BmGeneStart) %>% 
  mutate(BmGeneNum = row_number()) %>% ungroup()

maxGnBm.df <- gn.df %>% group_by(BmChr) %>% slice(which.max(BmGeneNum)) %>% 
  select(BmChr, BmGeneNum) %>% arrange(BmChr) %>% ungroup() %>% rename(BmGeneNumMax = BmGeneNum)

maxGnPn.df <- gn.df %>% group_by(PnChr) %>% slice(which.max(PnGeneNum)) %>% 
  select(PnChr, PnGeneNum) %>% arrange(PnChr) %>% ungroup() %>% rename(PnGeneNumMax = PnGeneNum)

condensed.df <- gn.df %>% left_join(maxGnBm.df, by = 'BmChr') %>% 
  left_join(maxGnPn.df, by = 'PnChr') %>% 
  mutate(BmGene.dfe = ifelse(BmGeneNum < BmGeneNumMax/2,
                             BmGeneNum, (BmGeneNumMax - BmGeneNum))) %>% 
  mutate(PnGene.dfe = ifelse(PnGeneNum < PnGeneNumMax/2,
                             PnGeneNum, (PnGeneNumMax - PnGeneNum))) %>% 
  select(BmChr, BmGeneStart, BmGeneNum, BmGeneNumMax, BmGene.dfe,
         PnChr, PnGeneStart, PnGeneNum, PnGeneNumMax, PnGene.dfe)

## Observed distribution of 
p1 <- ggplot(condensed.df) +
  geom_bin2d(aes(BmGene.dfe/BmGeneNumMax, PnGene.dfe/PnGeneNumMax), bins = 25) +
  viridis::scale_fill_viridis(limits = c(0,60),
                              guide_colorbar(raster = T)) +
  labs(x = "% Distance from end, Gene Order- Pnapi", y = '% Distance from end, Gene Order - Bmori') +
  guides(fill=guide_legend(title="Genes"))  +
  scale_x_continuous(breaks = c(0,.1,.2,.3,.4,.5),
                     labels = c(0,10,20,30,40,50),
                     expand = c(0,0)) +
  scale_y_continuous(breaks = c(0,.1,.2,.3,.4,.5),
                     labels = c(0,10,20,30,40,50),
                     expand = c(0,0)) +
  coord_fixed() +
  theme_minimal() +
  theme(plot.background = element_rect(fill = 'green', colour = 'red'))

q1 <- ggplot_build(p1)

result2 <- result
ggplot(result2) +
  geom_contour(aes(x,y,z = dir_p2))

## Estimated distribution plot
avgEpo10 <- Epo10 %>% group_by(x,y) %>% summarise(avgEst = mean(Estimated))

ggplot(avgEpo10) +
  geom_tile(aes(x,y, fill = avgEst)) +
  viridis::scale_fill_viridis()


### Remove last 3 smallest chromosomes and analyse distribution
largeChr <- BmoriPnapiSynt %>% filter(PnChr != 'Pnapi_chr_25' &
                                        PnChr != 'Pnapi_chr_24' &
                                        PnChr != 'Pnapi_chr_23' &
                                        PnChr != 'Pnapi_chr_22')

LcEpo10 <- data_frame()
for (i in 1:100){
  test <- reordChrom(largeChr) %>%
    arrange(fuseChr, fuseGeneNum)
  p <-test %>% filter() %>% ggplot() +
    geom_bin2d(aes(BmGene.dfe/BmGeneNumMax, fuseGene.dfe/fuseGeneNumMax), color = "black", bins = 10) +
    viridis::scale_fill_viridis() +
    theme_minimal() +
    labs(x = "% Distance from end - fuse", y = '% Distance from end - Bmori') +
    guides(fill=guide_legend(title="Genes")) +
    coord_fixed()
  q <- ggplot_build(p)
  LcEpo10 <- q$data %>% data.frame() %>% select(xbin, ybin, count) %>% 
    rename(Estimated = count, x = xbin, y = ybin) %>% bind_rows(LcEpo10)
}

## Combine estimated bins and average for Estimated probability plot
avgLcEpo10 <- LcEpo10 %>% group_by(x,y) %>% summarise(avgEst = mean(Estimated))

