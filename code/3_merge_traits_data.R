library(tidyverse)
library(ggpmisc)
library(ggridges)


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


checkm <- read.delim("data/CheckM_all.tsv", sep = "\t")
checkm <- checkm %>%
  select(BinId, Completeness, Contamination)

checkm$BinId <- gsub(checkm$BinId, pattern = ".fa", replacement = "")  
checkm <- checkm %>%
  rename(MAG = BinId)

taxon_data <- read.delim("data/GTDB_for_NCBItax.tsv", sep = " ", header = F)
colnames(taxon_data) <- c("taxonomy", "MAG")
taxon_data <- taxon_data %>%
  mutate(domain = str_extract(string = taxonomy, "(d_\\w+)"),
         phylum = str_extract(string = taxonomy, "(p_\\w+)"),
         family = str_extract(string = taxonomy, "(f_\\w+)")
  )

mag_strategies <- read.csv("data/mag_strategies.csv")
ENC <- read.delim("data/ENC_all_output.txt", sep = "\t", header = F) %>%
  select(-V6)
colnames(ENC) <- c("MAG", "ENC", "riboENC", "deltaENC", "ribo_cts")

AFE <- read.csv("data/AFE_4thwedge.csv", header = T)
AFE$MAG <- AFE$bin

size_and_gc <- read.csv("data/gc_and_size.csv", header = F, col.names = c("MAG", "gc", "bp", 
                                                                          "all_A", "all_C",
                                                                          "all_G", "all_T"))

nuc_freq_ribo <- read.csv("data/ribo_nucl_freq.csv", header = T)

nuc_freq_ribo_sum <- nuc_freq_ribo %>%
  rename(MAG = fasta)%>%
  group_by(MAG)%>%
  summarise(ribo_A_mean = weighted.mean(A, w = total_nuc),
            ribo_C_mean = weighted.mean(C, w = total_nuc),
            ribo_G_mean = weighted.mean(G, w = total_nuc),
            ribo_T_mean = weighted.mean(T, w = total_nuc),
            ribo_gene_length = mean(total_nuc))%>%
  mutate(ribo_GC = (ribo_G_mean+ribo_C_mean)/
           (ribo_A_mean+ribo_C_mean+ribo_G_mean+ribo_T_mean),
         ribo_GC_skew = (ribo_G_mean-ribo_C_mean)/(ribo_G_mean+ribo_C_mean),
         ribo_AT_skew = (ribo_A_mean-ribo_T_mean)/(ribo_A_mean+ribo_T_mean))

aa_ribo <- read_csv( "data/sum_ribo_aa.csv")

degen <- read_csv( "data/ribo_degenerate.csv")%>%
  select(-...1)

mag_strategies_w_traits <- full_join(mag_strategies, ENC) %>%
  full_join(., size_and_gc) %>%
  full_join(., taxon_data) %>%
  left_join(., nuc_freq_ribo_sum) %>%
  left_join(., aa_ribo)%>%
  left_join(., degen)

mag_strategies_w_traits$strategy_broad <- factor(mag_strategies_w_traits$strategy_broad, 
                                                 levels = c("Early Responders", "Middle Responders", 
                                                            "Late Responders", "Sensitive"))





summary(lm(riboENC ~ strategy_broad, mag_strategies_w_traits %>%
             drop_na(strategy_broad)))

AFE <- full_join(AFE, mag_strategies_w_traits) 

AFE <- left_join(AFE, checkm)

AFE <- AFE %>%
  mutate(adjusted_size = bp/(Completeness/100),
         ape.boot.CI = ape.boot.median - ape.boot.CI.L)

ns <- read.csv("data/ns_substitutions.csv") %>%
  select(-X)

mag_strategies_w_traits <- left_join(mag_strategies_w_traits, ns)
AFE <- left_join(AFE, ns)
AFE <- AFE %>% unique()

AFE$Mbp <- AFE$adjusted_size/1000000

write_csv(mag_strategies_w_traits, file = "data/MAG_traits_and_strategies.csv")


write_csv(AFE, file = "data/full_traits_growth_df.csv")
