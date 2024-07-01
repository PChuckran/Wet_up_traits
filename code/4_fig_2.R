library(tidyverse)
library(multcompView)
library(patchwork)



theme_pete<- function() {
  theme_minimal() %+replace%
    theme(
      text=element_text(size=12, family="Arial"),
      #axis.title=element_text(family="Futura Medium"),
      #axis.text=element_text(family="Optima"),
      strip.background =element_blank(),
      strip.placement = "outside",
      panel.grid.major.y = element_line(colour = "lightgrey"),
      panel.grid.minor.y = element_line(colour = "lightgrey"),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank()
    )
}

full_df <- read.csv("data/full_traits_growth_df.csv")%>%
  unique()

rna_df <- full_df %>% select(MAG, strategy_broad, riboENC, nsGC_skew,
                   degen_GCskew, nsAT_skew, degen_ATskew)%>%
  drop_na(strategy_broad) %>%
  unique()

hawks <- c("#002244", "#69BE28", "#A8DDFF", "#A5ACAF")
whawks <- c("#002244", "#69BE28", "#A8DDFF", "#A5ACAF")

full_df$strategy_broad <- factor(full_df$strategy_broad, levels = c("Early Responders",
                                                                    "Middle Responders",
                                                                    "Late Responders",
                                                                    "Sensitive"))


mod1 <- aov(riboENC~strategy_broad, full_df %>%
  #filter(ENC < 50)%>%
  select(strategy_broad, riboENC)%>%
  drop_na()%>%
  unique())

summary(mod1)
TukeyHSD(mod1)

fig2a <- full_df %>%
  #filter(ENC < 50)%>%
  select(strategy_broad, riboENC)%>%
  drop_na()%>%
  unique()%>%
  ggplot(., aes(riboENC, fill = strategy_broad, color = strategy_broad))+
  geom_density(alpha = 0.9)+
  scale_fill_manual("Transcriptional\nResponse", values = whawks)+
  scale_color_manual("Transcriptional\nResponse", values = hawks)+
  xlab("Effective Number of Codons (ENC`)\nof Ribosomal Protein Genes")+
  theme_pete()

fig2a



generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$type=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$type) , ]
  return(Tukey.labels)
}

traits_gc_skew <- rna_df %>%
  pivot_longer(cols = c(degen_GCskew, nsGC_skew),
               names_to = "type",
               values_to = "GC_skew" )

name_Switchgc <- Vectorize(function(a) {
  switch(as.character(a),
         'degen_GCskew' = '4-fold degenerate',
         'nsGC_skew' = 'Non-synonymous')
})


traits_gc_skew <- traits_gc_skew%>%
  mutate(type2 = name_Switchgc(type))


traits_gc_skew <- traits_gc_skew %>%
  mutate(strat_and_type = paste(strategy_broad, type))


codon_deg_mod <- aov(GC_skew ~ strat_and_type, traits_gc_skew)
codon_deg_anova <- anova(codon_deg_mod)
codon_deg_anova
f2b_tukey <- TukeyHSD(codon_deg_mod,  "strat_and_type", conf.level=0.95)
f2b_tukey <- generate_label_df(f2b_tukey, "strat_and_type")
f2b_tukey <- f2b_tukey %>%
  rename(strat_and_type = type)

f2b_tukey <- traits_gc_skew %>%
  select(strat_and_type, type2, strategy_broad, type)%>%
  unique()%>%
  left_join(f2b_tukey, .)

f2b <- traits_gc_skew %>%
  ggplot(., aes(strategy_broad, GC_skew, fill = type2))+
  # geom_text(data = f2b_tukey, aes(x = strategy_broad, y = .5, label = Letters), 
  #           inherit.aes = FALSE, color = "darkblue", size = 7)+
  theme_pete()+
  geom_hline(yintercept = 0, color = "darkred", size = 1, linetype = "dashed")+
  geom_boxplot(alpha = 0.95, width = .5)+
  scale_fill_manual("Site", values = c("#AA336A", "#3384AA"))+
  geom_text(data = f2b_tukey, aes(x = strategy_broad, y = .3, label = Letters, color = type2), 
            inherit.aes = FALSE, size = 5, position=position_dodge(width=0.7))+
  scale_color_manual("Site", values = c("#AA336A", "#3384AA"))+
  ylab("GC-skew")+
  xlab("Transcriptional response")+
  theme(plot.margin  = margin(.25, .25, .25, .25, "cm"))+
  coord_flip()+
  theme(legend.position = "none")

f2b


traits_at_skew <- rna_df %>%
  pivot_longer(cols = c(degen_ATskew, nsAT_skew),
               names_to = "type",
               values_to = "AT_skew" )

name_Switchat <- Vectorize(function(a) {
  switch(as.character(a),
         'degen_ATskew' = '4-fold degenerate',
         'nsAT_skew' = 'Non-synonymous')
})

traits_at_skew <- traits_at_skew%>%
  mutate(type2 = name_Switchat(type))


traits_at_skew <- traits_at_skew %>%
  mutate(strat_and_type = paste(strategy_broad, type))

codon_deg_mod <- aov(AT_skew ~ strat_and_type, traits_at_skew)
codon_deg_anova <- anova(codon_deg_mod)
codon_deg_anova
f2c_tukey <- TukeyHSD(codon_deg_mod,  "strat_and_type", conf.level=0.95)
f2c_tukey <- generate_label_df(f2c_tukey, "strat_and_type")
f2c_tukey <- f2c_tukey %>%
  rename(strat_and_type = type)

f2c_tukey <- traits_at_skew %>%
  select(strat_and_type, type2, strategy_broad, type)%>%
  unique()%>%
  left_join(f2c_tukey, .)


f2c <- traits_at_skew %>%
  ggplot(., aes(strategy_broad, AT_skew, fill = type2))+
  theme_pete()+
  geom_hline(yintercept = 0, color = "darkred", size = 1, linetype = "dashed")+
  geom_boxplot(alpha = 0.95, width = .5)+
  scale_fill_manual("Site", values = c("#AA336A", "#3384AA"))+
  geom_text(data = f2c_tukey, aes(x = strategy_broad, y = .4, label = Letters, color = type2), 
            inherit.aes = FALSE, size = 5, position=position_dodge(width=0.7))+
  scale_color_manual("Site", values = c("#AA336A", "#3384AA"))+
  ylab("AT-skew")+
  xlab("")+
  theme(plot.margin  = margin(.25, .25, .25, .25, "cm"))+
  coord_flip()+
  theme(axis.text.y=element_blank())
f2c

traits_at_skew

layout <- "
AA
BC"

Fig2 <- fig2a+f2b+f2c+
  plot_layout(design = layout)+
  plot_annotation(tag_levels = "A")
Fig2
ggsave(plot = Fig2, filename = "figures/Fig_2.svg", height = 8, width = 8)
