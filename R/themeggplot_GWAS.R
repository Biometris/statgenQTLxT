my_theme_empty_GWAS <-  
  theme(plot.background = element_blank(),
        panel.grid.major= element_line(colour="white"),
        panel.grid.minor= element_blank(),
        plot.title = element_text(size=20, face="bold", vjust=2) ,
        axis.ticks = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_blank(),
        legend.position="none",
        panel.background=element_rect(fill = "gray"),
        axis.text.x = element_text(size = 0),
        axis.text.y = element_text(size=13 , color = "black"), 
        axis.title.x = element_text(size=20, color = "navyblue"), 
        axis.title.y = element_text(size=20, color = "navyblue"), 
        strip.background=element_rect(fill = "gray40"),
        strip.text=element_text(),
        strip.text.x=element_text(size=14),
        strip.text.y=element_text(size=0) )

