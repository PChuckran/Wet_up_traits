## Calculate 4-fold degeneracy for ribosomal protein genes

library(tidyverse)

codon_aa_key <- read_csv("data/codon_aa_key.csv")

codon_aa_key <- codon_aa_key %>%
  rename(codon = Codon)%>%
  select(codon, AA)

ribo_cts <- read_csv("data/ribo_codon_freq.csv")%>%
  rename(MAG= fasta)

ribo_cts_long <- ribo_cts %>%
  pivot_longer(cols = -c(ID, MAG), names_to = "codon", values_to = "ct")%>%
  left_join(., codon_aa_key)

ribo_cts_sum <- ribo_cts_long %>%
  group_by(MAG, AA, codon)%>%
  summarise(sum_cts = sum(ct, na.rm = T))

ribo_cts_sum<- ribo_cts_long %>%
  group_by(MAG, AA)%>%
  summarise(all_aa_cts = sum(ct, na.rm = T))%>%
  left_join(ribo_cts_sum, .)

ribo_cts_sum <- ribo_cts_sum %>%
  mutate(codon_freq = sum_cts/all_aa_cts)

degen_AA <- c("A", "G", "P", "T", "V")

degen_ribo <- ribo_cts_sum%>%
  filter(AA %in% degen_AA)

degen_ribo <- degen_ribo %>%
  mutate(deg_site = str_sub(codon, start = 3, end = 3))

degen_ribo_sum <- degen_ribo %>%
  group_by(MAG, deg_site)%>%
  summarise(degen_ct = sum(sum_cts))

degen_ribo_sum_wide <- degen_ribo_sum %>%
  pivot_wider(names_from = deg_site, values_from = degen_ct)

degen_ribo_sum_wide <- degen_ribo_sum_wide %>%
  mutate(degen_GC = (G+C)/(A+C+G+T),
         degen_GCskew = (G-C)/(G+C),
         degen_ATskew = (A-T)/(A+T))

degen_ribo_sum_wide <- degen_ribo_sum_wide %>%
  select(-c(A,C,G,T))

write.csv(degen_ribo_sum_wide, file = "data/ribo_degenerate.csv")
