library(tidyverse)
library(ggpmisc)
library(patchwork)

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

my.formula <- y ~ poly(x, 1, raw=TRUE)




sub_df <- full_df %>%
  select(ENC, gc, riboENC, ribo_GC, MAG)%>%
  unique()

my.formula <- y ~ poly(x, 1, raw=TRUE)

F3a<- sub_df %>%
  ggplot(., aes(gc, ENC))+
  geom_point()+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(rr.label, ..p.value.label.., sep = "~~~")),
               small.p = T, 
               label.x = "right",
               label.y = "top",
               parse = TRUE,
               size = 2.8)+
  geom_smooth(method = "lm",  formula = my.formula, se = F,
              color = "black", linetype = "twodash")+
  xlab("Genome GC-%")+
  ylab("Genome codon usage bias (ENC`)")+
  theme_pete()

F3b<-sub_df %>%
  ggplot(., aes(ribo_GC*100, riboENC))+
  geom_point()+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(rr.label, ..p.value.label.., sep = "~~~")),
               small.p = T, 
               label.x = "right",
               label.y = "top",
               parse = TRUE,
               size = 2.8)+
  geom_smooth(method = "lm",  formula = my.formula, se = F,
              color = "black", linetype = "twodash")+
  xlab("Ribosomal protein gene GC-%")+
  ylab("Ribosomal protein gene\ncodon usage bias (ENC`)")+
  theme_pete()

Fs3 <- F3a+F3b

ggsave("figures/Fig_S3.svg", plot = Fs3, width = 8, height = 4)

full_df %>%
  drop_na(ape.boot.median)%>%
  ggplot(., aes(ribo_GC, riboENC, color = ape.boot.median))+
  geom_point()
