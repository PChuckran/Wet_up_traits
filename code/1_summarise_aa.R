## Calculate the C:N and total cost (in ~P bonds) of the AA in the ribosomal proteins 

library(tidyverse)
library(ggpmisc)

ribo_aa <- read.csv("data/ribo_aa_freq_cleaned.csv")%>%
  rename(MAG = fasta)

sum_ribo_aa <- ribo_aa %>%
  select(-c(ID, c_to_n))%>%
  group_by(MAG)%>%
  summarise_all(funs(sum))

sum_ribo_aa <- sum_ribo_aa %>%
  mutate(AA_total = A+R+N+D+C+E+Q+G+H+I+L+K+M+F+P+S+T+W+Y+V,
         AA_total_cost = (A*11.7)+(R*27.3)+(N*14.7)+(D*12.7)+(C*24.7)+
         (E*15.3)+(Q*16.3)+(G*11.7)+(H*38.3)+(I*32.3)+(L*27.3)+(K*30.3)+
         (M*34.3)+(F*52)+(P*20.3)+(S*11.7)+(T*18.7)+(W*74.3)+(Y*50)+(V*23.3),
         ribo_AA_mean_cost = AA_total_cost/AA_total) %>%
  mutate(ribo_aa_cn = carbon/nitrogen)


sum_ribo_aa <- sum_ribo_aa %>% select(MAG, ribo_aa_cn, ribo_AA_mean_cost)

write_csv(sum_ribo_aa, "data/sum_ribo_aa.csv")
