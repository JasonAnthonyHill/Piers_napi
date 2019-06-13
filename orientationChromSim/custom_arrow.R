my_geom_segment<- function (mapping = NULL, data = NULL, stat = "identity", position = "identity", 
                            ..., arrow = NULL, lineend = "butt", na.rm = FALSE, show.legend = NA, 
                            inherit.aes = TRUE) 
{
  layer(data = data, mapping = mapping, stat = stat, 
        geom = MyGeomSegment,                                       ###### <- changed this!
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(arrow = arrow, lineend = lineend, na.rm = na.rm, 
                      ...))
}

MyGeomSegment <- ggproto("GeomSegment", Geom,
                         required_aes = c("x", "y", "xend", "yend"),
                         non_missing_aes = c("linetype", "size", "shape"),
                         default_aes = aes(colour = "black", size = 0.5, linetype = 1, alpha = NA),
                         
                         draw_panel = function(data, panel_params, coord, arrow = NULL,
                                               lineend = "butt", na.rm = FALSE) {
                           
                           data <- remove_missing(data, na.rm = na.rm,
                                                  c("x", "y", "xend", "yend", "linetype", "size", "shape"),
                                                  name = "geom_segment")
                           if (ggplot2:::empty(data)) return(zeroGrob())        #### added ggplot2:::
                           
                           if (coord$is_linear()) {
                             coord <- coord$transform(data, panel_params)
                             return(segmentsGrob(coord$x, coord$y, coord$xend, coord$yend,
                                                 default.units = "native",
                                                 gp = gpar(
                                                   col = alpha(coord$colour, coord$alpha),
                                                   fill = alpha(coord$colour, coord$alpha),
                                                   lwd = coord$size * .pt,
                                                   lty = coord$linetype,
                                                   lineend = lineend,
                                                   linejoin = 'mitre'     #### <- added this!
                                                 ),
                                                 arrow = arrow
                             ))
                           }
                           
                           data$group <- 1:nrow(data)
                           starts <- subset(data, select = c(-xend, -yend))
                           ends <- plyr::rename(subset(data, select = c(-x, -y)), c("xend" = "x", "yend" = "y"),
                                                warn_missing = FALSE)
                           
                           pieces <- rbind(starts, ends)
                           pieces <- pieces[order(pieces$group),]
                           
                           GeomPath$draw_panel(pieces, panel_params, coord, arrow = arrow,
                                               lineend = lineend)
                         },
                         
                         draw_key = draw_key_path
)