library(tidyverse)
library(ggpmisc)


full_df <- read.csv("data/full_traits_growth_df.csv")

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

full_df$timef <- paste(full_df$Time,"h")
full_df$timef <- factor(x = full_df$timef, levels = c("24 h", "48 h", "72 h", "168 h"), ordered = T)
my.formula <- y ~ poly(x, 1, raw=TRUE)

Fig_S6 <- full_df %>%
  drop_na(Time)%>%
  ggplot(., aes(ribo_AA_mean_cost, ape.boot.median, color = timef))+
  geom_point()+
  stat_poly_eq(formula = my.formula,
               aes(label = paste(rr.label, ..p.value.label.., sep = "~~~")),
               small.p = T,
               label.x = "right",
               label.y = "top",
               parse = TRUE,
               size = 2.8)+
  geom_smooth(method = "lm", alpha = 0.5, se = FALSE, size = .8)+
  scale_color_viridis_d("Time")+
  xlab("Ribosomal protein mean amino acid cost (~P bonds)")+
  ylab("Atom Fraction Excess (AFE)")+
  theme_pete()

Fig_S6
ggsave(filename = "figures/Fig_S6.svg", Fig_S6, width = 5, height = 4)


