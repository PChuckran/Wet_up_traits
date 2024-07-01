library(tidyverse)
library(ggpmisc)


full_df <- read.csv("data/full_traits_growth_df.csv")

theme_border<- function() {
  theme_linedraw() %+replace%
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



full_df$timef <- paste(full_df$Time,"h")
full_df$timef <- factor(x = full_df$timef, levels = c("24 h", "48 h", "72 h", "168 h"), ordered = T)
my.formula <- y ~ poly(x, 1, raw=TRUE)


FS5a <- full_df %>%
  drop_na(Time)%>%
  ggplot(., aes(degen_GCskew, ape.boot.median, color = timef))+
  theme_border()+
  geom_point()+
  stat_poly_eq(formula = my.formula,
               aes(label = paste(rr.label, ..p.value.label.., sep = "~~~")),
               small.p = T,
               label.x = "right",
               label.y = "top",
               parse = TRUE,
               size = 2.8)+
  geom_smooth(method = "lm", alpha = 0.5, se = FALSE, size = .8)+
  geom_vline(xintercept = 0, color = "darkred", size = 1, linetype = "dashed")+
  xlab("Synonymous GC-skew")+
  ylab("AFE")+
  scale_color_viridis_d("Time")+
  scale_fill_viridis_d("Time")

FS5b <- full_df %>%
  drop_na(Time)%>%
  ggplot(., aes(degen_ATskew, ape.boot.median, color = timef))+
  theme_border()+
  geom_point()+
  stat_poly_eq(formula = my.formula,
               aes(label = paste(rr.label, ..p.value.label.., sep = "~~~")),
               small.p = T,
               label.x = "right",
               label.y = "top",
               parse = TRUE,
               size = 2.8)+
  geom_smooth(method = "lm", alpha = 0.5, se = FALSE, size = .8)+
  geom_vline(xintercept = 0, color = "darkred", size = 1, linetype = "dashed")+
  ylab("")+
  xlab("Synonymous AT-skew")+
  scale_color_viridis_d("Time")+
  scale_fill_viridis_d("Time")

FS5c <- full_df %>%
  drop_na(Time)%>%
  ggplot(., aes(nsGC_skew, ape.boot.median, color = timef))+
  theme_border()+
  geom_point()+
  geom_smooth(method = "lm", formula = my.formula, alpha = 0.5, se = FALSE, size = .8)+
  stat_poly_eq(formula = my.formula,
               aes(label = paste(rr.label, ..p.value.label.., sep = "~~~")),
               small.p = T,
               label.x = "left",
               label.y = "top",
               parse = TRUE,
               size = 2.8)+
  geom_vline(xintercept = 0, color = "darkred", size = 1, linetype = "dashed")+
  xlab("Non-synonymous GC-skew")+
  ylab("AFE")+
  scale_color_viridis_d("Time")+
  scale_fill_viridis_d("Time")


FS5d <- full_df %>%
  drop_na(Time)%>%
  ggplot(., aes(nsAT_skew, ape.boot.median, color = timef))+
  theme_border()+
  geom_point()+
  geom_smooth(method = "lm", alpha = 0.5, se = FALSE, size = .8)+
  stat_poly_eq(formula = my.formula,
               aes(label = paste(rr.label, ..p.value.label.., sep = "~~~")),
               small.p = T,
               label.x = "left",
               label.y = "top",
               parse = TRUE,
               size = 2.8)+
  geom_vline(xintercept = 0, color = "darkred", size = 1, linetype = "dashed")+
  xlab("Non-synonymous AT-skew")+
  ylab("")+
  scale_color_viridis_d("Time")+
  scale_fill_viridis_d("Time")

summary(lm(ape.boot.median~ Mbp+riboENC+ribo_GC+degen_GCskew+Time, full_df))
summary(lm(ape.boot.median~ degen_ATskew+Time, full_df))

Fig_S5 <- FS5a + FS5b + FS5c + FS5d +
  plot_layout(guides = 'collect')+
  plot_annotation(tag_levels = "A")

Fig_S5

ggsave(filename = "figures/Fig_S5.svg", width = 9, height = 7.5)
