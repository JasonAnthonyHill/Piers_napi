library(ggplot2)
library(data.table)
library(cowplot)
library(gridExtra)
library(grid)
library(scales)
library(dplyr)
library(gtable)
library(ggrepel)
source("custom_arrow.R")
#library(plyr)
setwd("/mnt/griffin/jashil/GenManFig_2")

set.seed(2017)

scafChunks <- fread("ScaffChrom_id.txt")
colnames(scafChunks) <- c("Chr", "Scaffold", "ID", "Start","End")

# syntBlocks <- data.frame(fread("ColorBmSynt.txt"), stringsAsFactors = T)
# colnames(syntBlocks) <- c("BmChr", "PnChr", "PnStart", "PnEnd", "hexCol")
# syntBlocks$BmChr <- factor(syntBlocks$BmChr, levels = c(paste(rep("Bmori_chr_",28),seq(1:28),sep = "")))

directionalSyntBlocks <- fread("directionalSyntBlocks.csv")

#syntLabels <- read.table("syntLabels.tsv")

mpSpanCov <- fread("mergedMPcov.tsv")
# rollReps.100k <- fread("rollSum/rollReps100k.tsv", header = T)

intraSNPs <- read.table(file = "intraChrSNP_l1.txt")
colnames(intraSNPs) <- c("ChromA","CoordA","ChromB","CoordB")

pasiLink <- fread(file = "pasiLinkageCoords.txt", col.names = c("mkNum", "Chr", "ChrPos", "cM"))

gp <- fread(file = "gp.tsv")

custPal <- fread(file = "BmChrColors.txt", header = F)
mycolors <- custPal$V2
names(mycolors) <- custPal$V1

## Subset data for quick testing
# mpSpanCovAll <- mpSpanCov
# mpSpanCov <- mpSpanCov[Chr == "Chromosome_2"][1:500000]
# rollReps.100kAll <- rollReps.100k
# rollReps.100k <- rollReps.100k[chr == "Chromosome_2"][1:500000]

## Individual plots used for figure
scaling <- 2

for (i in 1:25){
  # i <- 3
  # Scaffolds of Pn genome
  k <- paste("Chromosome_",i,sep = "")
  ks <- paste("Pnapi_chr_",i,sep = "")
  
  #scale limits by chromosome size
  xmax <- max(gp[gp$Chr == k,]$pos)
  scaleX <- scale_x_continuous(labels = c(0, round(0.25*xmax/1e6), round(0.5*xmax/1e6), round(0.75*xmax/1e6), round(xmax/1e6)), 
                               limits = c(1,xmax), breaks = c(0, round(0.25*xmax, digits = -6) , 
                                                              round(0.5*xmax, digits = -6), 
                                                              round(0.75*xmax, digits = -6), xmax))
  
  scafs.plot <- ggplot(data=scafChunks[scafChunks$Chr == k,]) +
    geom_rect(aes(xmin=Start,
                  ymin=1,
                  xmax=End,
                  ymax=2,
                  fill=Scaffold),
              color="black") +
    geom_label(aes(x=(Start+End)/2, y=1.5, label=ID), size=4.5*scaling, fill = 'white')+
    theme(axis.text=element_text(size=12*scaling),
          axis.line.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=12*scaling,face="bold",angle=0,vjust=1),
          legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank()) +
    scale_x_continuous(labels = c(0, round(0.25*xmax/1e6), round(0.5*xmax/1e6), round(0.75*xmax/1e6), round(xmax/1e6)), 
                       limits = c(1,xmax), breaks = c(0, round(0.25*xmax, digits = -6) , 
                                                      round(0.5*xmax, digits = -6), 
                                                      round(0.75*xmax, digits = -6), xmax)) +
    viridis::scale_fill_viridis(discrete = TRUE) +
    labs(y = "E)")
  scafs.gtable <- ggplot_gtable(ggplot_build(scafs.plot))
  
  # Syntenic genes between Pn and Bm
  synt.plot <- directionalSyntBlocks %>% filter(PnChr == ks) %>% 
    ggplot() +
    my_geom_segment(aes(x = ifelse(direction == "+", PnStart, PnEnd),
                        xend = ifelse(direction == "+", PnEnd*.95, PnStart*1.05),
                        y = 0,
                        yend = 0,
                        color = BmChr),
                    size = 6 * scaling,
                    arrow = arrow(length = unit(10 * scaling, "points"))) +
    ggrepel::geom_label_repel(aes(x = ifelse(BmCoord1 < BmCoord2,
                                             PnStart+(PnEnd-PnStart)*.4,
                                             PnStart+(PnEnd-PnStart)*.6),
                                  y = 0, 
                                  label = BmID),
                              force = .01,
                              direction = "x",
                              fill = "skyblue",
                              size = 3 * scaling) +
    theme(axis.text=element_text(size=12 * scaling),
          axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=12*scaling,face="bold",angle=0,vjust=.5),
          legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank()) +
    scaleX +
    scale_fill_manual(values = mycolors) +
    scale_color_manual(values = mycolors) +
    labs(y="D)")
  synt.gtable <- ggplot_gtable(ggplot_build(synt.plot))
  
  # Mate pair span coverage
  mp.plot <- ggplot(data = mpSpanCov[mpSpanCov$Chr == k]) +
    geom_line(aes(x = pos, y = sumMP)) +
    scale_y_log10() +
    #    ggtitle("Mate Pair Spanning Coverage") +
    theme(axis.text=element_text(size=12*scaling),
          axis.line.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=12*scaling,face="bold",angle=0,vjust=1),
          legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(),
          plot.title = element_text(hjust = 1)) +
    scaleX +
    labs(y = "A)")
  mp.gtable <- ggplot_gtable(ggplot_build(mp.plot))
  
  ## WGS SNPs
  chromSNPS <- intraSNPs[intraSNPs$ChromA == k & intraSNPs$ChromB == k, 1:2]
  uniqSNPs <- chromSNPS[!duplicated(chromSNPS),]
  wgssnp.plot <- ggplot(data = uniqSNPs) +
    geom_vline(xintercept = uniqSNPs$CoordA, color = "blue") +
    scale_y_continuous(limits = c(-1,1)) +
    scaleX +     
    theme(axis.text=element_text(size=12*scaling),
          axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=12*scaling,face="bold",angle=0,vjust=1),
          #          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank()) +
    labs(y = "C)")
  
  wgssnp.gtable <- ggplot_gtable(ggplot_build(wgssnp.plot))
  
  ## Pasi SNPs
  pasiLink.plot <- ggplot(data = pasiLink[Chr == k,]) +
    geom_point(aes(x = ChrPos, y = cM)) +
    theme(axis.text=element_text(size=12*scaling),
          axis.line.x=element_blank(),axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=12*scaling,face="bold",angle=0,vjust=1),
          legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank()) +
    scaleX +
    labs(y = "B)")
  pasiLink.gtable <- ggplot_gtable(ggplot_build(pasiLink.plot))
  
  
  ## Figure formatting
  # title <- ggdraw() + draw_label(paste("Chromosome",i,"- Pieris napi"), fontface='bold')
  maxWidth = unit.pmax(mp.gtable$widths[2:3],
                       synt.gtable$widths[2:3],scafs.gtable$widths[2:3],
                       pasiLink.gtable$widths[2:3],wgssnp.gtable$widths[2:3])
  
  mp.gtable$widths[2:3] <- maxWidth
  synt.gtable$widths[2:3] <- maxWidth
  scafs.gtable$widths[2:3] <- maxWidth
  pasiLink.gtable$widths[2:3] <- maxWidth
  wgssnp.gtable$widths[2:3] <- maxWidth
  
  fn <- paste("napiChromosome_",i,".png",sep = "")
  png(fn, width = 600 * scaling, height = 200 * scaling)
  # item <- cowplot::plot_grid(title, mp.plot,pasiLink.plot,wgssnp.plot,synt.p,scafs.plot,
  #                      ncol=1,
  #                      align = "v",
  #                      rel_heights = c(2,2,.75,2,1))
  item <- grid.arrange(mp.gtable,pasiLink.gtable,wgssnp.gtable,synt.gtable,scafs.gtable,
                       top = textGrob(paste("Chromosome",i,"- Pieris napi"), gp=gpar(fontsize=12*scaling)),
                       ncol=1,
                       heights = c(2,2,.45,1,1.25))
  
  print(item)
  dev.off()
}

