### THEME GRAY ####

theme_gray_etal<- function(base_size = 12, bgcolor = NA) 
{
  half_line <- base_size/2
  theme(
    line = element_line(colour = "black", size = rel(1.5), 
                        linetype = 1, lineend = "butt"), 
    rect = element_rect(fill = NA, colour = "black",
                        size = 0.5, linetype = 1),
    text = element_text(face = "plain",
                        colour = "black", size = base_size,
                        lineheight = 0.9,  hjust = 0.5,
                        vjust = 0.5, angle = 0, 
                        margin = margin(), debug = FALSE), 
    
    axis.line = element_blank(), 
    axis.text = element_text(size = rel(1.5), colour = "grey10"),
    axis.text.x = element_text(margin = margin(t = half_line/2), 
                               vjust = 1), 
    axis.text.y = element_text(margin = margin(r = half_line/2),
                               hjust = 1),
    axis.ticks = element_line(colour = "black", size=1), 
    axis.ticks.length = unit(half_line*0.75, "pt"), 
    axis.title = element_text(size = rel(1.5), colour = "black"),
    axis.title.x = element_text(margin = margin(t = half_line*5,
                                                b = half_line)),
    axis.title.y = element_text(angle = 90, 
                                margin = margin(r = half_line*5,
                                                l = half_line)),
    
    legend.background = element_rect(colour = NA), 
    legend.key = element_rect(colour = NA),
    legend.key.size = unit(2, "lines"), 
    legend.key.height = NULL,
    legend.key.width = NULL, 
    legend.text = element_text(size = rel(1)),
    legend.text.align = NULL,
    legend.title = element_text(size = rel(1)), 
    legend.title.align = NULL, 
    legend.position = "right", 
    legend.direction = NULL,
    legend.justification = "center", 
    legend.box = NULL, 
    
    panel.background = element_rect(fill=bgcolor,colour = "black", size = 2), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.spacing = unit(half_line, "pt"), panel.margin.x = NULL, 
    panel.spacing.y = NULL, panel.ontop = FALSE, 
    
    #Facet Labels
    strip.background = element_blank(),
    strip.text = element_text(face="bold",colour = "black", size = rel(1.5)),
    strip.text.x = element_text(margin = margin(t = half_line,
                                                b = half_line)), 
    strip.text.y = element_text(angle = 0, 
                                margin = margin(l = half_line, 
                                                r = half_line)),
    strip.switch.pad.grid = unit(5, "lines"),
    strip.switch.pad.wrap = unit(5, "lines"), 
    
    
    plot.background = element_rect(colour = rgb(119,136,153, max = 255)), 
    plot.title = element_text(size = rel(1.5), 
                              margin = margin(b = half_line * 1.2)),
    plot.margin = margin(4*half_line, 4*half_line, 4*half_line, 4*half_line),
    complete = TRUE)
}

