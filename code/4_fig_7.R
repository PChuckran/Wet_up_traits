library(tidyverse)
library(patchwork)
library(ggsci)
library(ggpmisc)

full_df <- read.csv("data/full_traits_growth_df.csv")%>%
  unique()

full_df <- full_df %>%
  filter(Contamination < 10)

full_df$strategy_broad <- factor(full_df$strategy_broad, levels = c("Early Responders",
                                                                    "Middle Responders",
                                                                    "Late Responders",
                                                                    "Sensitive"))

theme_pete<- function() {
  theme_bw() %+replace%
    theme(
      text=element_text(size=10, family="Arial"),
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

FS7a <- full_df %>%
  select(riboENC, MAG, phylum)%>%
  unique()%>%
  drop_na()%>%
  ggplot(., aes(phylum, riboENC))+
  geom_jitter(color = "darkgrey")+
  geom_boxplot(fill = "grey", alpha = 0.7)+
  theme_pete()+
  xlab("Phylum")+
  ylab("Ribosomal protein gene\ncodon bias (ENC`)")+
  theme(plot.margin = margin(0,1,0,0, "cm"))+
  coord_flip()


short <- c("#71A6F2", "#A6A77D")

full_df <- full_df %>%
  mutate(phylum = substr(phylum, 4, nchar(phylum)))


my.formula <- y ~ poly(x, 1, raw=TRUE)

FS7b <- full_df %>%
  #drop_na(Time)%>%
  select(MAG, Mbp, phylum, ape.boot.median, Time)%>%
  drop_na()%>%
  unique()%>%
  #filter(riboENC < 50)%>%
  #filter(Mbp < 10)%>%
  filter(phylum == "Actinobacteriota" | phylum == "Proteobacteria")%>%
  ggplot(., aes(Mbp, ape.boot.median, color = phylum))+
  geom_point()+
  #scale_color_brewer(palette = "Dark2")+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(rr.label, ..p.value.label.., sep = "~~~")),
               small.p = T, 
               label.x = "right",
               label.y = "top",
               parse = TRUE,
               size = 2.8)+
  theme_pete()+
  scale_color_manual(values = short)+
  scale_x_continuous(trans='log2')+
  xlab("Genome size (log2(Mbp))")+
  ylab("Atom fraction excess (AFE)")+
  geom_smooth(method = "lm",  formula = my.formula, se = F)+
  theme(legend.position = "bottom")+
  facet_wrap(.~Time)

FS7c <- full_df %>%
  #drop_na(Time)%>%
  select(MAG, degen_GCskew, phylum, strategy_broad)%>%
  #filter(riboENC < 50)%>%
  drop_na()%>%
  unique()%>%
  filter(phylum == "Actinobacteriota" | phylum == "Proteobacteria")%>%
  ggplot(., aes(strategy_broad, degen_GCskew, fill = phylum))+
  geom_boxplot()+
  #scale_color_brewer(palette = "Dark2")+
  theme_pete()+
  coord_flip()+
  scale_fill_manual(values = short)+
  ylab("Synonymous GC-skew")+
  xlab("Transcriptional response")+
  theme(legend.position = "none")+
  facet_wrap(.~phylum, scales = "free_x")

FS7d <- full_df %>%
  #drop_na(Time)%>%
  select(MAG, degen_ATskew, phylum, strategy_broad)%>%
  #filter(riboENC < 50)%>%
  drop_na()%>%
  unique()%>%
  filter(phylum == "Actinobacteriota" | phylum == "Proteobacteria")%>%
  ggplot(., aes(strategy_broad, degen_ATskew, fill = phylum))+
  geom_boxplot()+
  #scale_color_brewer(palette = "Dark2")+
  theme_pete()+
  coord_flip()+
  scale_fill_manual(values = short)+
  ylab("Synonymous AT-skew")+
  xlab("Transcriptional response")+
  theme(legend.position = "none")+
  facet_wrap(.~phylum, scales = "free_x")

FS7e <- full_df %>% 
  drop_na(Time)%>%
  unique()%>%
  filter(family == "f__Sphingomonadaceae" | family == "f__Burkholderiaceae")%>%
  ggplot(., aes(ribo_GC, ape.boot.median, color = family))+
  geom_point(data = full_df %>% drop_na(Time), 
             aes(ribo_GC, ape.boot.median), 
             color = "darkgrey", alpha = 0.8)+
  geom_point()+
  #scale_color_brewer(palette = "Dark2")+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(rr.label, ..p.value.label.., sep = "~~~")),
               small.p = T, 
               label.x = "right",
               label.y = "top",
               parse = TRUE,
               size = 2.8)+
  theme_pete()+
  scale_color_manual(values = c("red", "darkgreen"))+
  xlab("Ribosomal protein GC content")+
  ylab("Atom fraction excess (AFE)")+
  geom_smooth(method = "lm",  formula = my.formula, se = F)+
  theme(plot.margin = margin(0,0,0,0, "cm"))+
  theme(legend.position = "bottom")+
  
  facet_wrap(.~Time, scales = "free_x")






  
layout <- "
AC
AD
BE
"
Fig_S7 <- FS7a+FS7b+FS7c+FS7d+FS7e+plot_layout(design = layout)+
  plot_annotation(tag_levels = "A")

ggsave(filename = "figures/Fig_S7.svg", Fig_S7, height = 11.5, width = 9.5)
