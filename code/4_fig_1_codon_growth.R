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
      text=element_text(size=14, family="Arial"),
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

# f1a <- full_df %>%
#   #filter(Time == 24)%>%
#   drop_na(Time)%>%
#   ggplot(., aes(riboENC, ape.boot.median/Time, color = as.factor(Time)))+
#   geom_point(alpha = 0.9, size = 3)+
#   theme_pete()+
#   stat_poly_eq(formula = my.formula, 
#                aes(label = paste(rr.label, ..p.value.label.., sep = "~~~")),
#                small.p = T, 
#                label.x = "right",
#                label.y = "top",
#                parse = TRUE,
#                size = 2.8)+
#   #scale_color_viridis_d()+
#   scale_color_manual(values = As_pal)+
#   #geom_smooth(method = "lm", se = FALSE, formula = my.formula)+
#   geom_smooth(method = "lm", color = "black", formula = my.formula)+
#   facet_grid(Time~.)

f1a <- full_df %>%
  #filter(Time == 24)%>%
  filter(domain == "d__Bacteria")%>%
  #filter(ENC < 50)%>%
  drop_na(Time)%>%
  ggplot(., aes(riboENC, ape.boot.median, color = adjusted_size))+
  geom_point(alpha = 0.9, size = 3)+
  theme_pete()+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(rr.label, ..p.value.label.., sep = "~~~")),
               small.p = T, 
               label.x = "right",
               label.y = "top",
               parse = TRUE,
               size = 2.8)+
  scale_color_viridis_c()+
  #scale_color_manual(values = As_pal)+
  geom_smooth(method = "lm",  formula = my.formula, se = F)+
  facet_grid(Time~.)

f1a
my.formula <- y ~ poly(x, 1, raw=TRUE)
  
  full_df %>%
    filter(Time == "168 h")%>%
    #filter(domain == "d__Bacteria")%>%
    #filter(ENC < 50)%>%
    drop_na(Time)%>%
    ggplot(., aes(adjusted_size, ape.boot.median))+
    geom_point(color = "darkred", alpha = 0.9, size = 3)+
    theme_pete()+
    ylab("Atom Fraction Excess (AFE)")+
    xlab("Estimated Genome Size (Mbp)")+
    stat_poly_eq(formula = my.formula, 
                 aes(label = paste(rr.label, ..p.value.label.., sep = "~~~")),
                 small.p = T, 
                 label.x = "right",
                 label.y = "top",
                 parse = TRUE,
                 size = 2.8)+
    scale_color_viridis_c()+
    #scale_color_manual(values = As_pal)+
    geom_smooth(method = "lm",  formula = my.formula, se = F, color = "black")
  
  #ggsave("preliminary_figures/size_v_growth.png", width = 4, height = 3)
  
  my.formula <- y ~ poly(x, 1, raw=TRUE)
  
full_df %>%
  unique()%>%
  #filter(Time == 24)%>%
  #filter(domain == "d__Bacteria")%>%
  #filter(ENC < 50)%>%
  drop_na(Time)%>%
  ggplot(., aes(riboENC, ape.boot.median, color = Time))+
  geom_smooth(method = "lm",  formula = my.formula, se = F)+
  geom_point(alpha = 0.9, size = 3)+
  theme_pete()+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(rr.label, ..p.value.label.., sep = "~~~")),
               small.p = T, 
               label.x = "right",
               label.y = "top",
               parse = TRUE,
               size = 2.8)+
  scale_color_viridis_d()+
  ylab("Atom Fraction Excess (AFE)")+
  xlab("Ribosomal Protein
       Codon Usage Bias (ENC`)")

#ggsave("preliminary_figures/ribo_v_growth.png", width = 5.5, height = 4)


full_df %>%
  #filter(Time == 24)%>%
  #filter(domain == "d__Bacteria")%>%
  #filter(ENC < 50)%>%
  drop_na(Time)%>%
  ggplot(., aes(adjusted_size, ape.boot.median, color = Time))+
  geom_smooth(method = "lm",  formula = my.formula, se = F)+
  geom_point(alpha = 0.9, size = 3)+
  theme_pete()+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(rr.label, ..p.value.label.., sep = "~~~")),
               small.p = T, 
               label.x = "right",
               label.y = "top",
               parse = TRUE,
               size = 2.8)+
  scale_color_viridis_d(direction = -1)+
  ylab("Atom Fraction Excess (AFE)")

  scale_y_continuous(trans='log10')


name_Switch <- Vectorize(function(a) {
    switch(as.character(a),
           'riboENC' = 'ENC` (ribosomal protein genes)',
           'adjusted_size' = 'Estimated genome size (Mbp)',
           'ribo_GC' = 'GC (ribosomal protein genes')
  })



for_fig <- full_df %>%
  drop_na(Time)%>%
  filter(domain == "d__Bacteria")%>%
  pivot_longer(cols = c(riboENC, adjusted_size, ribo_GC), 
               values_to = "val", names_to = "trait")%>%
  mutate(trait = name_Switch(trait))

f2a <- for_fig %>%
  #filter(Time == 24)%>%
  #filter(ENC < 50)%>%
  ggplot(., aes(val, ape.boot.median))+
  geom_point(shape = 21, alpha = 0.8, size = 3, fill = "pink")+
  theme_pete()+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(rr.label, ..p.value.label.., sep = "~~~")),
               small.p = T, 
               label.x = "right",
               label.y = "top",
               parse = TRUE,
               size = 2.8)+
  xlab("")+
  ylab("AFE")+
  #scale_color_manual(values = As_pal)+
  geom_smooth(method = "lm",  formula = my.formula, se = F, color = "black")+
  facet_grid(Time~trait, scale = "free_x", switch = "x")

full_df %>%
  drop_na(Time)%>%
  filter(domain == "d__Bacteria")%>%
  ggplot(., aes(riboENC, ape.boot.median))+
  geom_point(size = 2)+
  ylab("Atom fraction excess")+
  xlab("Ribosomal protein\ncodon bias (ENC`)")+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(rr.label, ..p.value.label.., sep = "~~~")),
               small.p = T, 
               label.x = "right",
               label.y = "top",
               parse = TRUE,
               size = 2.8)+
  geom_smooth(method = "lm", se = F)+
  scale_color_manual(values = As_pal)+
  theme_pete()+
  facet_wrap(~Time, ncol = 1)



full_df %>%
  drop_na(Time)%>%
  filter(domain == "d__Bacteria")%>%
  ggplot(., aes(ribo_GC*100, ape.boot.median, color = riboENC))+
  geom_point( size = 2)+
  geom_point()+
  ylab("Atom fraction excess")+
  xlab("Ribosomal protein GC-%")+
  scale_color_gradient("Ribosomal\nprotein\nENC`",
                      low = "#EFB21E",
                      high = "#115740",
                      space = "Lab",
                      na.value = "grey50",
                      guide = "colourbar",
                      aesthetics = "color"
  )+ 
  theme_pete()+
  facet_wrap(.~Time, ncol = 1)

nar <- c("plum2", "coral1", "grey", "darkgreen")
As_pal <- c("black", "#4cbb17", "#EFB21E", "#A2AAAD")

full_df %>%
  drop_na(Time)%>%
  filter(domain == "d__Bacteria")%>%
  ggplot(., aes(ribo_GC*100, ape.boot.median, color = Time))+
  geom_point(aes(size = riboENC), alpha = 0.8)+
  ylab("Atom fraction excess")+
  xlab("Ribosomal protein GC-%")+
  scale_color_manual(values = As_pal)+
  #scale_size(trans = 'reverse')+
  theme_pete()+
  facet_wrap(.~Time)



my.formula <- y ~ poly(x, 1, raw=TRUE)

full_df %>%
  drop_na(Time)%>%
  filter(domain == "d__Bacteria")%>%
  select(Mbp, ape.boot.median, Time)%>%
  #filter(Mbp < 10)%>%
  unique()%>%
  ggplot(., aes(Mbp, ape.boot.median))+
  geom_point( alpha = 0.8, size =2)+
  ylab("Atom fraction excess")+
  xlab("Genome size (bp)")+
  # scale_color_manual(values = As_pal)+
  #scale_size(trans = 'reverse')+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(rr.label, ..p.value.label.., sep = "~~~")),
               small.p = T, 
               label.x = "right",
               label.y = "top",
               parse = TRUE,
               size = 2.8)+
  geom_smooth(method = "lm",  formula = my.formula, se = F)+
  theme_pete()+
  facet_wrap(.~Time, ncol = 1)




full_df %>%
  drop_na(Time)%>%
  filter(domain == "d__Bacteria")%>%
  ggplot(., aes(Mbp, ape.boot.median, color = riboENC))+
  geom_point( size = 2)+
  geom_point()+
  ylab("Atom fraction excess")+
 # xlab("Genome size (bp)")+
  scale_color_gradient("Ribosomal\nprotein\nENC`",
                       low = "#EFB21E",
                       high = "#115740",
                       space = "Lab",
                       na.value = "grey50",
                       guide = "colourbar",
                       aesthetics = "color"
  )+ 
  theme_pete()+
  facet_wrap(.~Time)






for_fig %>%
  filter(Time == 168)%>%
  #filter(ENC < 50)%>%
  ggplot(., aes(val, ape.boot.median))+
  geom_point(shape = 21, alpha = 0.8, size = 3, fill = "pink")+
  theme_pete()+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(rr.label, ..p.value.label.., sep = "~~~")),
               small.p = T, 
               label.x = "right",
               label.y = "top",
               parse = TRUE,
               size = 2.8)+
  xlab("")+
  ylab("AFE")+
  #scale_color_manual(values = As_pal)+
  geom_smooth(method = "lm",  formula = my.formula, se = F, color = "black")+
  facet_grid(Time~trait, scale = "free_x", switch = "x")

As_pal <- c("black", "#003831", "#EFB21E", "#A2AAAD")


f2b<- full_df %>%
  unique()%>%
  filter(domain == "d__Bacteria")%>%
  drop_na(Time)%>%
  filter(Time == "48 h")%>%
  ggplot(., aes(ribo_GC, ape.boot.median, fill = riboENC))+
  geom_point(shape = 21, alpha = 0.95, size = 3)+
  theme_pete()+
  scale_fill_gradient("Ribosomal\nprotein\nENC`",
    low = "#EFB21E",
    high = "#115740",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill"
  )+
  ylab("Atom Fraction Excess (AFE)")+
  xlab("GC content of\nribosomal protein genes")
f2b


ggsave(filename = "preliminary_figures/Fig1a.png", f2a, width = 5.5, height = 9)
ggsave(filename = "preliminary_figures/Fig1b.png", f2b, width = 4.5, height = 3.5)

full_df %>%
  filter(domain == "d__Bacteria")%>%
  drop_na(Time)%>%
  #filter(Time == "48 h")%>%
  ggplot(., aes(riboENC, ape.boot.median, fill = ribo_GC))+
  geom_point(shape = 21, alpha = 0.95, size = 3)+
  theme_pete()+
  scale_fill_gradient(,
                      low = "#EFB21E",
                      high = "#115740",
                      space = "Lab",
                      na.value = "grey50",
                      guide = "colourbar",
                      aesthetics = "fill"
  )+
  ylab("AFE")+
  xlab("GC content of\nribosomal protein genes")


summary(lm(ape.boot.median ~ ribo_GC+adjusted_size+riboENC, full_df))

