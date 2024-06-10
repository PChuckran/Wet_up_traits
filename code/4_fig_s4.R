library(tidyverse)

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

ribo_sum <- ribo %>%
  group_by(strategy_broad, Time)%>%
  summarise(LFC = mean(mean_lfc, na.rm = T),
            LFC_se = sd(mean_lfc, na.rm = T)/sqrt(length(MAG)))
 
t0 <- data_frame(strategy_broad = c("Early Responders", "Late Responders",   
                                    "Middle Responders", "Sensitive"),
                 Time = c(0,0,0,0),
                 LFC = c(0,0,0,0),
                 LFC_se = c(0,0,0,0))

ribo_sum <- rbind(ribo_sum, t0)

hawks <- c("#002244", "#69BE28", "#A8DDFF", "#A5ACAF", "white")

Fig_s4 <- ribo_sum %>%
  drop_na(strategy_broad)%>%
  ggplot(., aes(Time, LFC, color = strategy_broad))+
  geom_point(size = 2)+
  geom_errorbar(aes(ymin = LFC-LFC_se, 
                    ymax = LFC+LFC_se))+
  geom_line()+
  scale_color_manual("", values = hawks)+
  ylab("Mean upregulation of ribosomal
       protein genes (log2-fold change)")+
  theme_pete()

ggsave("figures/Fig_S4.svg", Fig_s4, width = 6, height = 4)
