### Relating ribosomal gene transcription to traits and growth ###
library(tidyverse)


taxon_data <- read.delim("data/GTDB_for_NCBItax.tsv", sep = " ", header = F)
colnames(taxon_data) <- c("taxonomy", "MAG")
taxon_data <- taxon_data %>%
  mutate(phylum = str_extract(string = taxonomy, "(p_\\w+)"),
         domain = str_extract(string = taxonomy, "(d_\\w+)"),
         family = str_extract(string = taxonomy, "(f_\\w+)"))

annotations <- read_delim("data/4w_annotations.tsv")
annotations <- annotations %>% dplyr::rename(Geneid = ...1,
                                             MAG = fasta)

checkm <- read.delim("data/CheckM_all.tsv", sep = "\t")
checkm <- checkm %>%
  select(BinId, Completeness, Contamination)

checkm$BinId <- gsub(checkm$BinId, pattern = ".fa", replacement = "")  
checkm <- checkm %>%
  rename(MAG = BinId)

de_data <- read.csv("data/4w_de_cts_FC.csv")%>%
  select(-X)

de_data <- de_data %>%
  mutate(time = str_extract(trt, "^([0-9]+)"),
         map = str_extract(trt, "([0-9]+$)"))


de_data <- left_join(de_data, annotations)

to_keep <- checkm %>% 
  filter(Contamination < 10)%>%
  filter(Completeness > 70)

de_data_hq <- de_data %>%
  filter(MAG %in% to_keep$MAG)

# if MAGs have fewer than n expressed genes, filter them out
threshold = 100
mag_counts <- as.data.frame(table(de_data_hq$MAG))
mags_to_keep <- mag_counts %>% filter(Freq > threshold)

de_data_hq <- de_data_hq %>%
  filter(MAG %in% mags_to_keep$Var1)

de_data_hq <- de_data_hq %>%
  mutate(reg = case_when(padj < 0.1 & log2FoldChange > 0 ~ "up", #
                         padj < 0.1 & log2FoldChange < 0 ~ "down", 
                         padj > 0.1 ~ "none",
                         time == "0" ~ "initial"
  ))

de_data_hq$time <- as.numeric(de_data_hq$time)

kegg_paths<- read.delim("data/ko_paths.txt", header = F)
colnames(kegg_paths) <- c("path", "kegg_id")
kegg_paths$kegg_id <- str_extract(kegg_paths$kegg_id, pattern = "K[0-9]+")


kegg_modules<- read.delim("data/module_to_ko.txt", header = F)
colnames(kegg_modules) <- c("kegg_id", "module")
kegg_modules <- kegg_modules %>%
  mutate(kegg_id = str_extract(kegg_id, pattern = "K[0-9]+"),
         module = str_extract(module, pattern = "M[0-9]+"))
module_names<- read.delim("data/module_names.txt", header = F)
colnames(module_names) <- c( "module", "module_name")
kegg_modules <- left_join(kegg_modules, module_names)

c_sources <- read.csv("data/C_sources_curated.csv")

strategies <- read.csv("data/mag_strategies.csv")

strategies <- strategies %>% select(MAG, strategy, strategy_broad)

de_data_hq <- left_join(de_data_hq, strategies)


de_data_hq$strategy_broad <- factor(de_data_hq$strategy_broad, 
                                    levels = c("Early Responders", "Middle Responders", 
                                               "Late Responders", "Sensitive"))

kegg_paths %>%
  filter(path == "path:map03010")%>%
  left_join(., de_data_hq) %>%
  #drop_na(padj)%>%
  drop_na(strategy_broad)%>%
  group_by(time, strategy_broad)%>%
  summarise(mean_lfc = mean(log2FoldChange, na.rm = T),
            se_lfc = sd(log2FoldChange, na.rm = T)/sqrt(length(log2FoldChange)))%>%
  ggplot(., aes(time, mean_lfc, color = strategy_broad))+
  geom_point()+
  geom_errorbar(aes(ymin=mean_lfc-se_lfc, ymax=mean_lfc+se_lfc))+
  geom_line()+
  scale_color_brewer("", palette = "Set1")+
  theme_minimal()+
  ylab("Log-2 Fold Change")

ribo_transcription_mags <- kegg_paths %>%
  filter(path == "path:map03010")%>%
  left_join(., de_data_hq) %>%
  filter(map == 100)%>%
  filter(baseMean >5)%>%
  group_by(MAG, strategy_broad, time)%>%
  summarise(base = mean(baseMean, na.rm = T),
            mean_lfc = mean(log2FoldChange, na.rm = T),
            se_lfc = sd(log2FoldChange, na.rm = T)/sqrt(length(log2FoldChange)))

ribo_transcription_mags %>%
  ggplot(., aes(time, mean_lfc, color = strategy_broad))+
  geom_point()+
  geom_errorbar(aes(ymin=mean_lfc-se_lfc, ymax=mean_lfc+se_lfc))+
  #geom_line()+
  scale_color_brewer("", palette = "Set1")+
  theme_minimal()+
  ylab("Log-2 Fold Change")

write.csv(ribo_transcription_mags, file = "ribo_transcription_mags.csv")

full_df <- read.csv(file = "data/full_traits_growth_df.csv")

moo <- left_join(ribo_transcription_mags %>% rename(Time = time), full_df)

my.formula <- y ~ poly(x, 1, raw=TRUE)

exp_v_growth <- moo %>%unique()%>%
  filter(Time != 3)%>%
  ggplot(. , aes(mean_lfc, ape.boot.median, color = as.factor(Time)))+
  ylab("Growth rate (atom fraction excess)")+
  xlab("Ribosomal Protein Gene Expression\n(log2-fold change)")+
  geom_point()+
  stat_poly_eq(formula = my.formula,
               aes(label = paste(rr.label, ..p.value.label.., sep = "~~~")),
               small.p = T,
               label.x = "left",
               label.y = "top",
               parse = TRUE,
               size = 2.8)+
  scale_color_brewer("Time", palette = "Set2")+
  geom_smooth(method = "lm", se = F)+
  theme_bw()

ggsave(plot = exp_v_growth, "preliminary_figures/exp_v_growth.png", 
       width = 5, height = 4)

summary(lm(ape.boot.median~mean_lfc+as.factor(Time), moo))

