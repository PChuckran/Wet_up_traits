library(tidyverse)
library(ggpmisc)
library(sjPlot)

full_df <- read.csv('data/full_traits_growth_df.csv')%>%
  unique() 

theme_pete<- function() {
  theme_bw() %+replace%
    theme(
      text=element_text(size=12, family="Arial"),
      #axis.title=element_text(family="Futura Medium"),
      #axis.text=element_text(family="Optima"),
      strip.background =element_blank(),
      strip.placement = "outside",
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank()
    )
}

full_df <- full_df %>%
  filter(Contamination < 10)

mod <- lm(ape.boot.median~riboENC+ribo_GC+adjusted_size+Time, full_df)
mod2 <- lm(ape.boot.median~riboENC*ribo_GC+adjusted_size+Time, full_df)
mod3 <- lm(ape.boot.median~riboENC*ribo_GC+Time, full_df)


AIC(mod, mod2, mod3)

FS2 <- plot_model(mod, type = "pred", terms = c("ribo_GC", "riboENC"))+
  scale_color_viridis_d("Ribosomal Protein\nCodon Bias (ENC`)")+
  scale_fill_viridis_d("Ribosomal Protein\nCodon Bias (ENC`)")+
  theme_pete()+
  xlab("Ribosomal Protein GC content")+
  ylab("Atom Fraction Excess (AFE)")+
  ggtitle("")

ggsave("figures/Fig_S2.svg", width = 5, height = 4)

  
  
