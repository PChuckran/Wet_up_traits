library(tidyverse)
library(ggExtra)
library(patchwork)
library(ggpmisc)


full_df <- read.csv('data/full_traits_growth_df.csv')%>%
  unique() 

full_df <- full_df %>%
  filter(Contamination < 10)

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

As_pal <- c("black", "#003831", "#EFB21E", "#A2AAAD")

full_df$Time <- paste(full_df$Time, "h", sep = " ")

full_df$Time <- factor(full_df$Time, level= c("24 h", "48 h", "72 h", "168 h"))

my.formula <- y ~ poly(x, 1, raw=TRUE)


f1a <- full_df %>%
  #filter(domain == "d__Bacteria")%>%
  drop_na(Time)%>%
  ggplot(., aes(riboENC, ape.boot.median))+
  geom_point(size = 3, shape = 21, fill = "lightgrey", alpha = 0.6)+
  theme_pete()+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(rr.label, ..p.value.label.., sep = "~~~")),
               small.p = T, 
               label.x = "right",
               label.y = "top",
               parse = TRUE,
               size = 2.8)+
  #scale_color_manual(values = As_pal)+
  geom_smooth(method = "lm",  formula = my.formula, se = F,
              color = "black", linetype = "twodash")+
  ylab("Atom Fraction Excess (AFE)")+
  xlab("Ribosomal Protein
       Codon Usage Bias (ENC`)")+
  facet_grid(Time~.)
f1a


f1b <- full_df %>%
  drop_na(Time)%>%
  #filter(Mbp < 10)%>%
  ggplot(., aes(Mbp, ape.boot.median))+
  geom_point(size = 3, shape = 21, fill = "lightgrey", alpha = 0.6)+
  theme_pete()+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(rr.label, ..p.value.label.., sep = "~~~")),
               small.p = T, 
               label.x = "right",
               label.y = "top",
               parse = TRUE,
               size = 2.8)+
  #scale_color_manual(values = As_pal)+
  geom_smooth(method = "lm",  formula = my.formula, se = F,
              color = "black", linetype = "twodash")+
  ylab("")+
  xlab("Estimated Genome Size (Mbp)")+
  facet_grid(Time~.)
f1b

f1c <- full_df %>%
  drop_na(Time)%>%
  filter(domain == "d__Bacteria")%>%
  ggplot(., aes(ribo_GC*100, ape.boot.median, color = riboENC))+
  geom_point( size = 2)+
  geom_point()+
  ylab("")+
  xlab("Ribosomal protein GC-%")+
  scale_color_gradient("Ribosomal\nprotein\nENC`",
                       low = "#EFB21E",
                       high = "darkblue",
                       space = "Lab",
                       na.value = "grey50",
                       guide = "colourbar",
                       aesthetics = "color"
  )+ 
  theme_pete()+
  facet_grid(Time~.)


Fig_1 <- f1a+f1b+f1c+  plot_annotation(tag_levels = "A")
Fig_1

ggsave(plot = Fig_1, filename = "figures/Fig_1.svg", width = 9, height = 8)





