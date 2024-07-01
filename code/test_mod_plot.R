library(tidyverse)
library(sjPlot)



full_df <- read.csv("data/full_traits_growth_df.csv")%>%
  unique()

mod <- lm(ape.boot.median~riboENC+ribo_GC+adjusted_size+Time, full_df)
mod2 <- lm(ape.boot.median~riboENC*ribo_GC+adjusted_size+Time, full_df)
mod3 <- lm(ape.boot.median~riboENC*ribo_GC+Time, full_df)


AIC(mod, mod2, mod3)

plot_model(mod)

summary(mod)
#  
# ggplot(full_df, aes(riboENC, ape.boot.median, color = ribo_GC))+
#   geom_point()+

plot_model(mod, type = "pred", terms = c("ribo_GC", "riboENC"))
  scale_color_viridis_d("Ribosomal\nGC")+
  scale_fill_viridis_d("Ribosomal\nGC")+
  theme_classic()
  xlab("Ribosomal Codon Bias (ENC`)")+
  ylab("Atom Fraction Excess (AFE)")+
  ggtitle("")


full_df %>%
  ggplot(., aes(gc, ENC, color = domain))+
  geom_point()+
  geom_smooth(method = "lm", se = F)

full_df %>%
  ggplot(., aes(ribo_GC, riboENC, color = domain))+
  geom_point()+
  geom_smooth(method = "lm", se = F)

full_df %>%
  ggplot(., aes(gc, ape.boot.median, color = ENC))+
  geom_point()+
  geom_smooth(method = "lm", se = F)

full_df %>%
  #filter(phylum == "p__Proteobacteria")%>%
  drop_na(ape.boot.median)%>%
  ggplot(., aes(ribo_GC, riboENC, color = ape.boot.median))+
  geom_point(size = 3)+
  scale_color_viridis_c()+
  facet_grid(phylum~Time)

full_df %>%
  #filter(phylum == "p__Proteobacteria")%>%
  drop_na(ape.boot.median)%>%
  ggplot(., aes( degen_GC, ape.boot.median, color = riboENC))+
  geom_point(size = 3)+
  scale_color_viridis_c()
  facet_grid(phylum~Time)

full_df %>%
  #filter(phylum == "p__Proteobacteria")%>%
  drop_na(ape.boot.median)%>%
  ggplot(., aes( gc, ape.boot.median, color = riboENC))+
  geom_point(size = 3)+
  scale_color_viridis_c(direction = -1)+
  facet_grid(.~Time)

full_df %>%
  mutate(gc_ribo_bias = ((ribo_GC*100)-gc)/gc)%>%
  drop_na(ape.boot.median)%>%
  ggplot(., aes(ribo_GC, gc_ribo_bias))+
  geom_point()+
  geom_smooth(method = "lm", se = F)+
  facet_wrap(.~Time)

full_df %>%
  drop_na(strategy_broad)%>%
  ggplot(., aes(strategy_broad, ENC))+
  geom_boxplot()

for_paper <- full_df %>%
  filter(Time == 168) %>%
  drop_na(ape.boot.median) #%>%
  #filter(adjusted_size < 8000000)

for_paper %>%
  ggplot(., aes(adjusted_size, ape.boot.median))+
  geom_point()

mod_for_paper <- lm(ape.boot.median ~ adjusted_size, for_paper)

to_predict <- data.frame(adjusted_size = c(2000000, 8000000))
predict(mod_for_paper, to_predict)

summary(lm(deltaENC~strategy_broad, full_df %>%
  drop_na(strategy_broad)))

summary(lm(ENC~gc_ribo_bias, full_df %>%
  mutate(gc_ribo_bias = ((ribo_GC*100)-gc)/gc)))

summary(lm(ENC~gc, full_df))
  
