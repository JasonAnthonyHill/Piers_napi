## Combine 4 block plots into single figure

library(ggplot2)
library(cowplot)
library(dplyr)

Opo10.gtable <- ggplotGrob(Opo10.plot) # convert to gtable
permEopo.gtable <- ggplotGrob(permEopoPlot) # convert to gtable
permEnd.gtable <- ggplotGrob(permEndPlot)
permGO.gtable <- ggplotGrob(permGOPlot)
p <- plot_grid(permEnd.gtable, permGO.gtable, 
               Opo10.gtable, permEopo.gtable, align = c('h'), ncol = 2,
               labels = 'auto')
p

p <- plot_grid(Opo10.gtable, permEopo.gtable, align = c('h'), ncol = 1)
p

png(filename = 'BlockOrderContent_panelsCD_geneOrder_Comparison.png', width = 1200, height = 1200)
p
dev.off()

png(filename = 'BlockOrderFig_C_geneOrder.png', width = 400, height = 400)
Opo10.plot
dev.off()

png(filename = 'BlockOrderFig_D_geneOrderPval.png', width = 400, height = 400)
permEopoPlot
dev.off()
