library(tidyverse)
library(cowplot)
library(ggsci)
library(ggpmisc)



full_df <- read.csv("data/full_traits_growth_df.csv")%>%
  unique()

full_df <- full_df %>%
  filter(Contamination < 10)

theme_pete<- function() {
  theme_bw() %+replace%
    theme(
      text=element_text(size=11, family="Arial"),
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

full_df$strategy_broad <- factor(full_df$strategy_broad, levels = c("Early Responders",
                                                                    "Middle Responders",
                                                                    "Late Responders",
                                                                    "Sensitive"))
giants <- c( "#EFD19F", "#FD5A1E", "#27251F", "white")

full_df <- full_df %>%
  mutate(phylum = substr(phylum, 4, nchar(phylum)))



ct_by_strat_phy <- full_df %>%
  select(MAG, strategy_broad, phylum)%>%
  unique()%>%
  group_by(strategy_broad, phylum)%>%
  summarise(ct = length(MAG))

ct_by_strat_phy <- full_df %>%
  select(MAG, strategy_broad, phylum)%>%
  unique()%>%
  group_by(strategy_broad)%>%
  summarise(ct_strat = length(MAG))%>%
  left_join(ct_by_strat_phy, .)%>%
  mutate(perc = (ct/ct_strat),
         type = "Transcribing Taxa")%>%
  rename(category = strategy_broad)

full_df$Timef <- factor(paste(full_df$Time, "h", sep = " "), levels = c("24 h", "48 h", "72 h", "168 h"))

ct_by_strat_phy_growth <- full_df %>%
  select(MAG, Timef, phylum, ape.boot.median)%>%
  drop_na()%>%
  group_by(Timef, phylum)%>%
  summarise(ct = length(MAG))

ct_by_strat_phy_growth <- full_df %>%
  select(MAG, Timef, phylum, ape.boot.median)%>%
  drop_na()%>%
  unique()%>%
  group_by(Timef)%>%
  summarise(ct_strat = length(MAG))%>%
  left_join(ct_by_strat_phy_growth, .)%>%
  mutate(perc = (ct/ct_strat),
         type = "Growing Taxa")%>%
  rename(category = Timef)

ct_by_strat_phy_growth$category <- as.factor(ct_by_strat_phy_growth$category)

ct_all <- rbind(ct_by_strat_phy, ct_by_strat_phy_growth)

moab <- c("#C19082","#71A6F2","#331A13","#7B6A62",
          "#466CAC","#5C5E2F", "#A6A77D", "#694335","#EAA683")

short <- c("#71A6F2", "#A6A77D")


F4A <- ct_all %>%
  drop_na(category)%>%
  ggplot(., aes(category, perc, fill = phylum))+
  geom_bar(stat = "identity")+
  theme_minimal()+
  xlab("")+
  #scale_fill_simpsons()+
  scale_fill_manual(values = moab)+
  #scale_fill_brewer("Phylum", palette = "Set3")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())+
  coord_flip()+
  scale_x_discrete(limits=rev)+
  facet_wrap(.~type, scales = "free", ncol = 1)+
  theme(legend.position = "right")

F4A

skoo <- full_df %>%
  filter(!is.na(ape.boot.median) | !is.na(strategy_broad))

F4B <- full_df %>%
  filter(!is.na(ape.boot.median) | !is.na(strategy_broad))%>%
  select(riboENC, MAG, phylum)%>%
  unique()%>%
  drop_na()%>%
  ggplot(., aes(phylum, riboENC, color = phylum, fill = phylum))+
  geom_jitter(alpha = 0.5, size = 1)+
  geom_boxplot(alpha = 0.7)+
  theme_pete()+
  scale_color_manual(values = moab)+
  scale_fill_manual(values = moab)+
  xlab("Phylum")+
  ylab("Ribosomal protein gene\ncodon bias (ENC`)")+
  coord_flip()+
  theme(legend.position = "none")+
  scale_x_discrete(limits=rev)+
  theme(legend.position = "none",
        plot.margin = margin(1,1,1,1, "cm"))

full_df %>%
  filter(!is.na(ape.boot.median) | !is.na(strategy_broad))%>%
  select(riboENC, MAG, phylum)%>%
  unique()%>%
  drop_na()%>%
  ggplot(., aes(phylum, riboENC, color = phylum, fill = phylum))+
  geom_jitter(alpha = 0.5, size = 1)+
  geom_violin(alpha = 0.7)+
  theme_pete()+
  scale_color_manual(values = moab)+
  scale_fill_manual(values = moab)+
  xlab("Phylum")+
  ylab("Ribosomal protein gene\ncodon bias (ENC`)")+
  coord_flip()+
  theme(legend.position = "none")+
  scale_x_discrete(limits=rev)



F4B

F4C <- full_df %>%
  #drop_na(Time)%>%
  select(MAG, riboENC, phylum, strategy_broad)%>%
  #filter(riboENC < 50)%>%
  drop_na()%>%
  unique()%>%
  filter(phylum == "Actinobacteriota" | phylum == "Proteobacteria")%>%
  ggplot(., aes(strategy_broad, riboENC, fill = phylum))+
  geom_boxplot()+
  #scale_color_brewer(palette = "Dark2")+
  theme_pete()+
  coord_flip()+
  scale_fill_manual(values = short)+
  ylab("Ribosomal protein gene\ncodon bias (ENC`)")+
  xlab("Transcriptional response")+
  facet_wrap(.~phylum)+
  theme(legend.position = "none")+
  theme(legend.position = "none",
        plot.margin = margin(0,1,0,1, "cm"))



my.formula <- y ~ poly(x, 1, raw=TRUE)

F4D <- full_df %>%
  #drop_na(Time)%>%
  select(MAG, riboENC, phylum, ape.boot.median, Time)%>%
  drop_na()%>%
  unique()%>%
  filter(riboENC < 50)%>%
  filter(phylum == "Actinobacteriota" | phylum == "Proteobacteria")%>%
  ggplot(., aes(riboENC, ape.boot.median, color = phylum))+
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
  scale_x_continuous(trans='log10')+
  xlab("Ribosomal protein gene\ncodon bias (ENC`)")+
  ylab("Atom fraction excess (AFE)")+
  geom_smooth(method = "lm",  formula = my.formula, se = F)+
  facet_wrap(.~Time)+
  theme(legend.position = "none",
        plot.margin = margin(0,2,0,2, "cm"))
          
F4D

 
 # layout <- "
 # AB
 # AB
 # CD
 # "
# 
# layout <- "
# A#
# AC
# BD
# B#
# "
# 
# Fig_4 <- F4A+F4B+F4D+F4C+
#    plot_layout(design = layout)+
#    plot_annotation(tag_levels = "a")
# Fig_4
# 
# Fig_4 <- plot_grid(F4A, F4B, F4D, F4C, labels = "B", rel_heights = c(2,1))
# fbcd <- plot_grid(F4B, F4D, F4C, labels = c("B", "C", "D"), ncol = 1, rel_heights = c(3,2,2))
# Fig_4 <- plot_grid(F4A, fbcd, ncol = 2, labels = c("A", ""), rel_widths = c(3,2))
# Fig_4

Fig_4 <- plot_grid(F4A, F4B, F4D, F4C, labels = c("A", "B", "C", "D"), ncol = 2, rel_heights = c(3,2,2))

ggsave(filename = "figures/Fig_4.png", width = 12, height = 11, Fig_4)



full_df %>%
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
  ylab("Ribosomal Protein Gene\nCodon Bias (ENC`)")+
  xlab("Transcriptional response")+
  facet_wrap(.~phylum, scales = "free_x")

full_df %>%
  #drop_na(Time)%>%
  select(MAG, ribo_AA_mean_cost, phylum, strategy_broad)%>%
  #filter(riboENC < 50)%>%
  drop_na()%>%
  unique()%>%
  filter(phylum == "Actinobacteriota" | phylum == "Proteobacteria")%>%
  ggplot(., aes(strategy_broad, ribo_AA_mean_cost, fill = phylum))+
  geom_boxplot()+
  #scale_color_brewer(palette = "Dark2")+
  theme_pete()+
  coord_flip()+
  scale_fill_manual(values = short)+
  ylab("Ribosomal Protein Gene\nCodon Bias (ENC`)")+
  xlab("Transcriptional response")+
  facet_wrap(.~phylum, scales = "free_x")




full_df %>%
  #drop_na(Time)%>%
  select(MAG, Mbp, phylum, ape.boot.median, Time)%>%
  drop_na()%>%
  unique()%>%
  filter(Mbp < 10)%>%
  #filter(riboENC < 50)%>%
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
  scale_x_continuous(trans='log10')+
  xlab("Genome Size (Mbp)")+
  ylab("Atom fraction excess (AFE)")+
  geom_smooth(method = "lm",  formula = my.formula, se = F)+
  facet_wrap(.~Time)

testy_df <- full_df %>%
  #drop_na(Time)%>%
  unique()%>%
  drop_na(ape.boot.median)%>%
  group_by(family, Time)%>%
  summarise(ct = length(MAG)) %>%
  left_join(full_df, .)

testy_df %>%
  #filter(ct >2)%>%
  #drop_na(Time)%>%
  unique()%>%
  #filter(riboENC < 50)%>%
  #filter( phylum == "Proteobacteria")%>%
  filter(family != "f__Sphingomonadaceae")%>%
  ggplot(., aes( ribo_GC, ape.boot.median, color = phylum))+
  geom_point()+
  #scale_color_brewer(palette = "Dark2")+
  # stat_poly_eq(formula = my.formula, 
  #              aes(label = paste(rr.label, ..p.value.label.., sep = "~~~")),
  #              small.p = T, 
  #              label.x = "right",
  #              label.y = "top",
  #              parse = TRUE,
  #              size = 2.8)+
  geom_smooth(method = "lm", se = F)+
  theme_pete()+
  facet_wrap(.~Time)

full_df %>%
  #drop_na(Time)%>%
  unique()%>%
  drop_na(Time)%>%
  #filter(Mbp < 10)%>%
  #filter(riboENC < 45)%>%
  filter(family == "f__Sphingomonadaceae")%>%
  ggplot(., aes(riboENC, ribo_GC, color = ape.boot.median))+
  geom_point()+
  #scale_color_brewer(palette = "Dark2")+
  # stat_poly_eq(formula = my.formula, 
  #              aes(label = paste(rr.label, ..p.value.label.., sep = "~~~")),
  #              small.p = T, 
  #              label.x = "right",
  #              label.y = "top",
  #              parse = TRUE,
  #              size = 2.8)+
  theme_pete()+
  #scale_color_manual(values = short)+
  #geom_smooth(method = "lm",  formula = my.formula, se = F)+
  facet_wrap(.~Time)

full_df %>%
  #drop_na(Time)%>%
  select(MAG, Mbp, phylum, ape.boot.median, ribo_GC, Time, riboENC)%>%
  drop_na()%>%
  unique()%>%
  #filter(Mbp < 10)%>%
  filter(riboENC < 50)%>%
  filter(phylum == "Actinobacteriota")%>%
  ggplot(., aes(ribo_GC, ape.boot.median, color = riboENC))+
  geom_point()+
  #scale_color_brewer(palette = "Dark2")+
  # stat_poly_eq(formula = my.formula, 
  #              aes(label = paste(rr.label, ..p.value.label.., sep = "~~~")),
  #              small.p = T, 
  #              label.x = "right",
  #              label.y = "top",
  #              parse = TRUE,
  #              size = 2.8)+
  theme_pete()+
  #scale_color_manual(values = short)+
  ylab("Atom fraction excess (AFE)")+
  geom_smooth(method = "lm",  formula = my.formula, se = F)+
  facet_wrap(.~Time)

summarized_df <- testy_df %>%
  filter(ct > 2)%>%
  #drop_na(Time)%>%
  select(MAG, Mbp, family, ape.boot.median, ribo_GC, Time, riboENC)%>%
  drop_na()%>%
  unique()%>%
  group_by(family, Time)%>%
  summarise_all(mean, na.rm = T)

summarized_df %>%
  ggplot(., aes(, ape.boot.median))+
  geom_point()+
  facet_wrap(.~Time)

full_df %>%
  #drop_na(Time)%>%
  select(MAG, riboENC, phylum, ribo_GC)%>%
  drop_na()%>%
  unique()%>%
  #filter(Mbp < 10)%>%
  #filter(riboENC < 50)%>%
  filter(phylum == "Actinobacteriota" | phylum == "Proteobacteria")%>%
  ggplot(., aes(riboENC, ribo_GC, color = phylum))+
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
  geom_smooth(method = "lm",  formula = my.formula, se = F)


temp <- full_df %>%
  #drop_na(Time)%>%
  select(riboENC, MAG, Mbp, phylum, ape.boot.median, ribo_GC, Time)%>%
  drop_na()%>%
  unique()%>%
  #filter(Mbp < 10)%>%
  filter(Time == 24)%>%
  #filter(riboENC < 50)%>%
  #filter(phylum == "Actinobacteriota")
  filter(phylum == "Proteobacteria")

summary(lm(ape.boot.median~riboENC+ribo_GC+Time, temp))


full_df %>%
  #drop_na(Time)%>%
  select(MAG, degen_GCskew, phylum, strategy_broad)%>%
  drop_na()%>%
  unique()%>%
  ggplot(., aes(strategy_broad, degen_GCskew, color = phylum))+
  geom_boxplot()+
  #scale_color_brewer(palette = "Dark2")+
  theme_pete()+
  coord_flip()+
  ylab("Codon Bias")+
  xlab("Transcriptional response")+
  facet_wrap(.~phylum, nrow = 1)

my.formula <- y ~ poly(x, 1, raw=TRUE)

summary(aov(riboENC ~ strategy_broad, full_df %>%
              #drop_na(Time)%>%
              select(MAG, riboENC, phylum, strategy_broad)%>%
              drop_na()%>%
              unique()%>%
              filter(phylum == "Proteobacteria")))

summary(aov(riboENC ~ strategy_broad, full_df %>%
              #drop_na(Time)%>%
              select(MAG, riboENC, phylum, strategy_broad)%>%
              drop_na()%>%
              unique()%>%
              filter(phylum == "Actinobacteriota")))

summary(aov(ribo_GC_skew ~ strategy_broad, full_df %>%
              #drop_na(Time)%>%
              select(MAG, ribo_GC_skew, phylum, strategy_broad)%>%
              drop_na()%>%
              unique()%>%
              filter(phylum == "Proteobacteria")))

summary(aov(ribo_GC_skew ~ strategy_broad, full_df %>%
              #drop_na(Time)%>%
              select(MAG, ribo_GC_skew, phylum, strategy_broad)%>%
              drop_na()%>%
              unique()%>%
              filter(phylum == "Actinobacteriota")))


