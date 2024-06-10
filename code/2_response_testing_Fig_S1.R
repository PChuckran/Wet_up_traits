library(tidyverse)

taxon_data <- read.delim("data/GTDB_for_NCBItax.tsv", sep = " ", header = F)
colnames(taxon_data) <- c("taxonomy", "MAG")
taxon_data <- taxon_data %>%
  mutate(phylum = str_extract(string = taxonomy, "(p_\\w+)"))

checkm <- read.delim("data/CheckM_all.tsv", sep = "\t")
checkm <- checkm %>%
  select(BinId, Completeness, Contamination)

checkm$BinId <- gsub(checkm$BinId, pattern = ".fa", replacement = "")  
checkm <- checkm %>%
  rename(MAG = BinId)


annotations <- read_delim("data/4w_annotations.tsv")
annotations <- annotations %>% dplyr::rename(Geneid = ...1,
                                             MAG = fasta)

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
  filter(MAG %in% to_keep$MAG) %>%
  drop_na(kegg_id)%>%
  drop_na(padj)

# if MAGs have fewer than n expressed genes, filter them out
threshold = 100
mag_counts <- as.data.frame(table(de_data_hq$MAG))
mags_to_keep <- mag_counts %>% filter(Freq > threshold)

kegg_de <- de_data_hq %>%
  filter(MAG %in% mags_to_keep$Var1)

kegg_de <- kegg_de %>%
  mutate(reg = case_when(padj < 0.1 & log2FoldChange > 0 ~ "up", #
                         padj < 0.1 & log2FoldChange < 0 ~ "down", 
                         padj > 0.1 ~ "none",
                         time == "0" ~ "initial"
  ))

mag_reg_totals <- kegg_de %>%
  group_by(MAG, time, reg) %>%
  summarise(genes = n())

mag_reg_totals <- mag_reg_totals %>%
  ungroup()%>%
  group_by(MAG, time)%>%
  summarise(total_genes = sum(genes))%>%
  left_join(mag_reg_totals, .)%>%
  mutate(rel_genes = genes/total_genes)

#mag_reg_totals_for_graph <- mag_reg_totals %>%
#  mutate(rel_genes = ifelse(reg == "down", -1*rel_genes, rel_genes))

mag_reg_totals <- mag_reg_totals %>%
  left_join(., taxon_data)

mag_reg_totals <- mag_reg_totals %>%
  mutate(family = str_extract(string = taxonomy, "(f_\\w+)"))

#mag_reg_totals$reg <- factor(mag_reg_totals$reg, levels = c("up", "none", "down"))

mag_reg_totals$time <- as.numeric(mag_reg_totals$time)
mag_reg_totals$timef <- factor(mag_reg_totals$time, 
                              levels = c("0", "3", "24", "48", "72", "168"))

exp_by_taxa <- mag_reg_totals %>%
  filter(reg != "none")%>%
  ggplot(., aes(timef, rel_genes*100, color = reg, group = paste(MAG, reg)))+
  geom_point()+
  geom_line()+
  ylab("Proportion of expressed genes (%)")+
  scale_color_brewer("Direction", palette = "Set1")+
  facet_wrap(paste(phylum, family, sep="\n")~.)+
  theme_linedraw(base_size = 10)+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())

exp_by_taxa

mag_reg_totals_wide <- mag_reg_totals %>%
  select(-genes)%>%
  pivot_wider(., names_from = reg, values_from = rel_genes)

mag_reg_totals_wide[is.na(mag_reg_totals_wide)] <- 0

mag_reg_totals_wide <- mag_reg_totals_wide %>%
  mutate(differential = up - down)

#mag_reg_totals_wide <-

mag_reg_totals_wide <- mag_reg_totals_wide %>%
  group_by(MAG) %>%
  filter(differential == max(differential))%>%
  rename(maxtime = time)%>%
  ungroup() %>%
  select(MAG, maxtime)%>%
  ungroup()%>%
  group_by(MAG)%>%
  # Is there more than one max? if so, set to NA
  summarise(max_time = mean(as.numeric(maxtime)),
            n_max = length(maxtime)) %>%
  mutate(max_time = ifelse(n_max > 1, "NA", max_time)) %>%
  select(MAG, max_time) %>%
  left_join(mag_reg_totals_wide, .)


mag_reg_diff <- mag_reg_totals_wide %>%
  select(MAG, time, differential) %>%
  pivot_wider(names_prefix = "time_", names_from = time, values_from = differential)

mag_reg_diff <- mag_reg_totals_wide %>%
  ungroup()%>%
  select(MAG, max_time) %>%
  unique()%>%
  left_join(mag_reg_diff, .)

mag_strategies <- mag_reg_diff %>%
  mutate(strategy_broad = case_when(
      max_time == 0 ~ "Sensitive",
      max_time == 3 ~ "Early Responders",
      max_time == 24 | max_time == 48 | max_time == 72  ~ "Middle Responders",
      max_time == 168 ~ "Late Responders")
  ) 



EO_prop_upreg <- mag_reg_totals_wide %>%
  left_join(., mag_strategies %>% select(MAG, strategy_broad)) %>%
  filter(strategy_broad == "Early Responders") %>%
  ggplot(., aes(timef, differential, group = MAG))+
  geom_point()+
  geom_line()+
  ylab("Net proportion of genes upregulated")+
  xlab("Time (h)")+
  facet_wrap(paste(phylum, family)~., scales = "free_y")+
  theme_linedraw(base_size = 10)+
  geom_hline(color = "red", yintercept = 0, linetype = "dashed")+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())



all_diff_reg <- mag_reg_totals_wide %>%
  left_join(., mag_strategies %>% select(MAG, strategy_broad)) %>%
  group_by(phylum, family, time) %>%
  summarise(differential = mean(differential, na.rm = T))%>%
  ggplot(., aes(time, differential))+
  geom_point()+
  geom_line()+
  ylab("Net proportion of genes upregulated")+
  xlab("Time (h)")+
  facet_wrap(paste(phylum, family)~., scales = "free_y")+
  theme_linedraw(base_size = 10)+
  scale_color_brewer(palette = "Set1")+
  geom_hline(color = "red", yintercept = 0, linetype = "dashed")+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())
all_diff_reg
#ggsave(filename = "preliminary_figures/egu/all_diff_reg.png", all_diff_reg, width = 18, height = 9)


mag_reg_totals_wide %>%
  left_join(., mag_strategies %>% select(MAG, strategy_broad)) %>%
  ggplot(., aes(time, differential, group = MAG, color = strategy_broad))+
  geom_point()+
  geom_line()+
  ylab("Net proportion of genes upregulated")+
  xlab("Time (h)")+
  facet_wrap(paste(phylum, family)~., scales = "free_y")+
  theme_linedraw(base_size = 10)+
  scale_color_brewer(palette = "Set1")+
  geom_hline(color = "red", yintercept = 0, linetype = "dashed")+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())


mag_reg_totals_wide %>%
  left_join(., mag_strategies %>% select(MAG, strategy_broad)) %>%
  ggplot(., aes(timef, differential, group = MAG, color = strategy_broad))+
  geom_point()+
  geom_line()+
  ylab("Net proportion of genes upregulated")+
  xlab("Time (h)")+
  facet_wrap(paste(phylum, family)~., scales = "free_y")+
  theme_linedraw(base_size = 10)+
  geom_hline(color = "red", yintercept = 0, linetype = "dashed")+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())


mag_strategies$strategy_broad <- factor(mag_strategies$strategy_broad, 
                                        levels = c("Early Responders", 
                                                   "Middle Responders",
                                                   "Late Responders", 
                                                   "Sensitive"))

strategies_overview <- mag_reg_totals_wide %>%
  left_join(., mag_strategies %>% select(MAG, strategy_broad)) %>%
  select(strategy_broad, time, differential) %>%
  ungroup()%>%
  group_by(strategy_broad, time)%>%
  summarise_all(mean)%>%
  drop_na(strategy_broad)%>%
  ggplot(., aes(time, differential, color = strategy_broad))+
  geom_point(size = 2)+
  geom_line()+
  ylab("Net proportion of genes upregulated")+
  xlab("Time (h)")+
  facet_wrap(strategy_broad~.)+
  theme_linedraw(base_size = 10)+
  scale_color_brewer(palette = "Set1")+
  geom_hline(color = "red", yintercept = 0, linetype = "dashed")+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none")

strategies_overview
ggsave(strategies_overview, filename = "figures/Fig_S1.png", 
       width = 6, height = 4.5)


  
write.csv(mag_strategies, file = "data/mag_strategies.csv")

mag_reg_totals_wide <- mag_reg_totals_wide %>%
  left_join(., mag_strategies %>% select(MAG, strategy_broad))

taxa_summaries <- mag_reg_totals_wide%>%
  ungroup()%>%
  select(MAG, phylum, strategy_broad)%>%
  unique(.)%>%
  group_by(phylum, strategy_broad)%>%
  summarise(count = length(MAG))

taxa_summaries %>%
  drop_na()%>%
  ggplot(., aes(count, phylum, fill = strategy_broad))+
  geom_bar(stat = "identity")+
  scale_fill_brewer(palette = "Dark2")

taxa_summaries <- taxa_summaries %>%
  drop_na()%>%
  group_by(strategy_broad)%>%
  summarise(total_strat = sum(count))%>%
  left_join(taxa_summaries, . )%>%
  mutate(perc_strat = count/total_strat)

taxa_summaries <- taxa_summaries %>%
  drop_na()%>%
  ungroup()%>%
  group_by(phylum)%>%
  summarise(total_phy = sum(count))%>%
  left_join(taxa_summaries, . )%>%
  mutate(perc_phy = count/total_phy)


taxa_summaries %>%
  drop_na()%>%
  ggplot(., aes(perc_strat, strategy_broad, fill = phylum))+
  geom_bar(stat = "identity")+
  theme_minimal()+
  ylab("")+
  xlab("")+
  
  scale_fill_brewer(palette = "Set1")

taxa_summaries %>%
  drop_na()%>%
  ggplot(., aes(perc_phy, phylum, fill = strategy_broad))+
  geom_bar(stat = "identity")+
  xlab("")+
  scale_fill_brewer(palette = "Dark2")

### IGNORE ANYTHING BELOW ###


mag_reg_totals %>%
  filter(reg != "none")%>%
  ggplot(., aes(time, genes, color = reg, group = paste(MAG, reg)))+
  geom_point()+
  geom_line()+
  #ylab("Proportion of expressed genes (%)")+
  scale_color_brewer("Direction", palette = "Set1")+
  facet_wrap(paste(phylum, family, sep="\n")~.)+
  theme_linedraw(base_size = 10)+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())


#ggsave(plot = exp_by_taxa, filename = "preliminary_figures/reg_by_taxa.pdf")
mag_reg_totals$time <- as.numeric(mag_reg_totals$time)
mag_reg_totals <- mag_reg_totals %>%
  filter(reg == "up")%>%
  group_by(MAG)%>%
  filter(rel_genes == max(rel_genes))%>%
  rename(maxtime = time) %>%
  ungroup()%>%
  select(MAG, maxtime)%>%
  group_by(MAG)%>%
  summarise(maxtime = mean(maxtime, na.rm = T))%>%
  left_join(mag_reg_totals, .)

mag_reg_totals <- mag_reg_totals %>%
  mutate(response_type = case_when(maxtime < 24 ~ "early", #
                                   maxtime > 71 ~ "late", 
                                   TRUE ~ "mid"
  ))

write_csv(mag_reg_totals, file = "data/mag_de_totals.csv")


  


