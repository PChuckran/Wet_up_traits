library(tidyverse)

codon_sites <- read.csv(file = "data/codon_degeneracy.csv")

codon_sites <- codon_sites %>%
  mutate(p1 = substr(degeneracy, 1, 1),
         p2 = substr(degeneracy, 2, 2),
         p3 = substr(degeneracy, 3, 3),
         c1 = substr(codon, 1, 1),
         c2 = substr(codon, 2, 2),
         c3 = substr(codon, 3, 3))

codon_sites <- codon_sites %>%
  mutate(c1 = ifelse(p1 != 1, NA, c1),
         c2 = ifelse(p2 != 1, NA, c2),
         c3 = ifelse(p3 != 1, NA, c3))


codon_sites_summary <- codon_sites %>%
  select(codon, c1, c2, c3) %>%
  pivot_longer(cols = c(c1, c2, c3)) 

codon_sites_summary_wide <- codon_sites_summary%>%
  drop_na()%>%
  select(-name)%>%
  table(.)%>%
  as.data.frame(.)%>%
  pivot_wider(names_from = value, values_from = Freq, names_prefix = "ns_")



ribo_cts <- read_csv("data/ribo_codon_freq.csv")%>%
  rename(MAG= fasta)  

ribo_cts <- ribo_cts %>%
  pivot_longer(cols = -c(MAG, ID), names_to = "codon", values_to = "count")

ribo_cts <- left_join(ribo_cts, codon_sites_summary_wide)

ribo_cts_summary <- ribo_cts %>%
  group_by(MAG)%>%
  summarise(nsA = sum(ns_A*count, na.rm = T),
            nsC = sum(ns_C*count, na.rm = T),
            nsG = sum(ns_G*count, na.rm = T),
            nsT = sum(ns_T*count, na.rm = T))%>%
  mutate(nstotal = nsA+nsC+nsG+nsT,
         nsA = nsA/nstotal,
         nsC = nsC/nstotal,
         nsG = nsG/nstotal,
         nsT = nsT/nstotal)

ribo_cts_summary <- ribo_cts_summary %>%
  mutate(nsGC = nsC+nsG,
         nsGC_skew = (nsG-nsC)/(nsG+nsC),
         nsAT_skew = (nsA-nsT)/(nsA+nsT)
         )

write.csv(ribo_cts_summary, file = "data/ns_substitutions.csv")

ribo_cts_ID_summary <-ribo_cts %>%
  group_by(ID)%>%
  summarise(nsA = sum(ns_A*count, na.rm = T),
            nsC = sum(ns_C*count, na.rm = T),
            nsG = sum(ns_G*count, na.rm = T),
            nsT = sum(ns_T*count, na.rm = T))%>%
  mutate(nstotal = nsA+nsC+nsG+nsT,
         nsA = nsA/nstotal,
         nsC = nsC/nstotal,
         nsG = nsG/nstotal,
         nsT = nsT/nstotal)

ribo_cts_ID_summary <- ribo_cts_ID_summary %>%
  mutate(nsGC = nsC+nsG,
         nsGC_skew = (nsG-nsC)/(nsG+nsC),
         nsAT_skew = (nsA-nsT)/(nsA+nsT)
  )

write.csv(ribo_cts_ID_summary, file = "data/ns_substitutions_by_gene.csv")

##this is weird##
ribo_cts_ID_summary %>%
  ggplot(., aes(nsGC, nsAT_skew))+
  geom_point()+
  geom_smooth(method = "lm")

ribo_cts_summary %>%
  ggplot(., aes(nsGC, nsAT_skew))+
  geom_point()

summary(lm(nsAT_skew~nsGC, ribo_cts_ID_summary))

test_skews <- codon_sites_summary_wide %>%
  mutate(at_skew = (ns_A-ns_T)/(ns_A+ns_T),
         gc_skew = (ns_G-ns_C)/(ns_G+ns_C),
         gc = (ns_G+ns_C)+(ns_G+ns_C+ns_A+ns_T))

## maybe this is it? 
## come back to this
test_skews %>%
  ggplot(., aes(gc, at_skew))+
  geom_point()+
  geom_smooth(method = "lm")


summary(lm(at_skew~gc, test_skews))
