library(tidyverse)



full_df <- read.csv("data/full_traits_growth_df.csv")%>%
  unique()

mod <- lm(ape.boot.median~ribo_GC+riboENC, full_df)

plot_model(mod)

full_df %>%
  filter(riboENC < 52)%>%
  filter(Time == 72)%>%
  ggplot(., aes(ribo_GC, ape.boot.median))+
  geom_point()+
  geom_smooth(method = "lm", se = F)
  facet_wrap(.~Time)
  
  full_df %>%
    filter(riboENC < 52)%>%
    filter(Time == 72)%>%
    ggplot(., aes(ribo_GC, riboENC))+
    geom_point()+
    geom_smooth(method = "lm", se = F)
  facet_wrap(.~Time)

hawks <- c("#002244", "#69BE28", "#A5ACAF", "#A8DDFF")
whawks <- c("#002244", "#69BE28", "#A5ACAF", "#A8DDFF")

full_df$strategy_broad <- factor(full_df$strategy_broad, levels = c("Early Responders",
                                                                    "Middle Responders",
                                                                    "Late Responders",
                                                                    "Sensitive"))
wgiants <- c( "#27251F", "#EFD19F", "white", "#FD5A1E")
bgiants <- c("#27251F", "#EFD19F", "black", "#FD5A1E")


fig2a <- full_df %>%
  select(strategy_broad, riboENC)%>%
  drop_na()%>%
  unique()%>%
  ggplot(., aes(riboENC, fill = strategy_broad, color = strategy_broad))+
  geom_density(alpha = 0.9)+
  scale_fill_manual("Transcriptional Response\nTiming", values = whawks)+
  scale_color_manual("Transcriptional Response\nTiming", values = hawks)+
  xlab("Effective Number of Codons (ENC`)")+
  theme_minimal()

fig2a

ggsave(plot = fig2a, filename = "preliminary_figures/Fig2a.png", width = 5, height = 3)


summarized_fd <- full_df %>%
  select(strategy_broad, Time, ape.boot.median)%>%
  drop_na()%>%
  unique()%>%
  group_by(strategy_broad, Time)%>%
  summarise(mean_ape = mean(ape.boot.median, na.rm = T),
            se_ape = sd(ape.boot.median, na.rm = T)/sqrt(n()))

summarized_fd %>%
  ggplot(., aes(Time, mean_ape, color = strategy_broad, fill = strategy_broad))+
  geom_point(shape = 21, size = 3)+
  geom_errorbar(aes(ymin=mean_ape-se_ape, ymax=mean_ape+se_ape, width = 1))+
  geom_line()+
  scale_fill_manual("Transcriptional Response\nTiming", values = whawks)+
  scale_color_manual("Transcriptional Response\nTiming", values = hawks)

full_df %>%
  select(strategy_broad, riboENC, degen_GCskew)%>%
  drop_na()%>%
  unique()%>%
  ggplot(., aes(riboENC, degen_GCskew, color = strategy_broad))+
  geom_point()+
  geom_smooth(method = "lm", color = "black")+
  theme_minimal()


summary(lm(riboENC~degen_ATskew+strategy_broad, full_df))

