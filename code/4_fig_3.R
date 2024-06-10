library(tidyverse)
library(ggpmisc)
library(ggsci)


full_df <- read.csv("data/full_traits_growth_df.csv")%>%
  unique()

ribo <- read.csv("data/ribo_transcription_mags.csv") %>%
  rename(Time = time)%>%
  select(-X)

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

full_df <- left_join(full_df, ribo)

full_df$Time <- as.factor(full_df$Time)

full_df$strategy_broad <- factor(full_df$strategy_broad, levels = 
                                   c("Early Responders",
                                     "Middle Responders",
                                     "Late Responders",
                                     "Sensitive"))

my.formula <- y ~ poly(x, 1, raw=TRUE)

time_pal <- c( "gold", "#FD5A1E", "cyan4", "darkblue" )


          
F3a <- full_df %>%
  drop_na(Time)%>%
  ggplot(., aes(mean_lfc, ape.boot.median, color = Time))+
  geom_point(size = 3, alpha = 0.9)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(rr.label, ..p.value.label.., sep = "~~~")),
               small.p = T, 
               label.x = "left",
               label.y = "top",
               parse = TRUE,
               size = 2.8)+
  ylab("Atom fraction excess (AFE)")+
  xlab(expression("Upregulation of ribosomal protein genes (lo" * g[2] * "-fold change)"))+
  scale_color_manual(values = time_pal)+
  geom_smooth(method = "lm", se = F)+
  theme_pete()

hawks <- c("#002244", "#69BE28", "#A8DDFF", "#A5ACAF", "white")


F3b <- full_df %>%
  drop_na(Time)%>%
  drop_na(strategy_broad)%>%
  ggplot(., aes(Time, ape.boot.median, fill = strategy_broad))+
  ylab("Atom fraction excess (AFE)")+
  xlab("Time after rewetting (h)")+
  scale_fill_manual("", values = hawks)+
  geom_boxplot( alpha = 0.7)+
  theme_pete()
  #theme(legend.position = "none")
F3b

summary(aov(ape.boot.median~strategy_broad, full_df %>%
              drop_na(Time)%>%
              drop_na(strategy_broad)%>%
              filter(Time == "24")
)
)

summary(aov(ape.boot.median~strategy_broad, full_df %>%
              drop_na(Time)%>%
              drop_na(strategy_broad)%>%
              filter(Time == "48")
)
)

summary(aov(ape.boot.median~strategy_broad, full_df %>%
              drop_na(Time)%>%
              drop_na(strategy_broad)%>%
              filter(Time == "72")
)
)

summary(aov(ape.boot.median~strategy_broad, full_df %>%
              drop_na(Time)%>%
              drop_na(strategy_broad)%>%
              filter(Time == "168")
)
)


summarized_df <- full_df %>%
  drop_na(Time)%>%
  drop_na(strategy_broad)%>%
  group_by(strategy_broad, Time)%>%
  summarise(ct = length(MAG))

summarized_df <-full_df %>%
  ungroup()%>%
  drop_na(Time)%>%
  drop_na(strategy_broad)%>%
  group_by( Time)%>%
  summarise(ct_total = length(MAG))%>%
  left_join(summarized_df, .)%>%
  mutate(perc = ct/ct_total)

F3c <- summarized_df %>%
  ggplot(., aes(Time, perc, fill = strategy_broad))+
  geom_bar(stat = "identity")+
  xlab("Time after rewetting (h)")+
  ylab("")+
  scale_fill_manual("", values = hawks)+
  theme_pete()

layout <- "
AB
#C
"

Fig_3 <- F3a+F3b+F3c+
  plot_layout(design = layout)+
  plot_annotation(tag_levels = "A")

Fig_3

ggsave(filename = "figures/Fig_3.svg", Fig_3, width = 10, height = 7)
