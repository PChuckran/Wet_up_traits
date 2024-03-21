library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(vegan)
library(tibble)
library(tidyverse)


theme_simple<- function() {
  theme_bw() %+replace%
    theme(
      text=element_text(size=10, family="Lucida Grande"),
      strip.background = element_rect(fill = "white"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank()
    )
}

# Get list of annotated geneids
annotations <- read_delim("data/4w_annotations.tsv")

annotations <- annotations %>% dplyr::rename(Geneid = ...1,
                                             MAG = fasta)
annotations <- annotations %>%
  mutate(is_annotated = case_when(
    is.na(kegg_id) & is.na(peptidase_id) & is.na(cazy_hits) &
      is.na(vogdb_categories) & is.na(pfam_hits)  ~ "no",
    TRUE ~ "yes"
  ))

# filter by annotated geneIds
annotated_geneIDs <- annotations %>% 
  filter(is_annotated == "yes") %>%
  select(Geneid)


# Load all data
all_data <- read.delim("data/old/all_cts_all_times.txt", sep = "\t", header = T)

rownames(all_data) <- all_data$Geneid

quality_cts_annotated <- all_data %>%
  filter(Geneid %in% annotated_geneIDs$Geneid | str_detect(Geneid, "tRNA"))

quality_cts_annotated <- quality_cts_annotated %>% select(starts_with("X"))






quality_cts_annotated[is.na(quality_cts_annotated)] <- 0


input_data <- as.matrix(quality_cts_annotated)

rownames(input_data) <- ref_annotations$Geneid

coldata <- data.frame(file_name = colnames(input_data)) %>%
  mutate(time = str_match(file_name, "12C18O_([0-9]+)H")[,2],
         map = str_match(file_name, "12C18O_[0-9]+H_([0-9]+)")[,2],
         sample = str_match(file_name, "(12C18O_[0-9]+H_[0-9]+_P[0-9]+_MT)")[,2])

coldata$time <- as.factor(coldata$time)
coldata$map <- as.factor(coldata$map)

setequal(coldata$file_name, colnames(input_data))
setdiff(coldata$file_name, colnames(input_data))


dds <- DESeqDataSetFromMatrix(countData = input_data,
                              colData = coldata,
                              design = ~ map+time+map:time)

dds$map <- relevel(dds$map, ref = "100")
dds$time <- relevel(dds$time, ref = "0")

keep <- rowSums(counts(dds)) > 100
dds <- dds[keep,]

dds <- DESeq(dds)
resultsNames(dds)

summary(results(dds, name = "time_3_vs_0"))



time_3<- as.data.frame(results(dds, name = "time_3_vs_0", test = "Wald"))
time_3 <- tibble::rownames_to_column(time_3, "Geneid")
time_3$trt <- "3_100"
time_24<- as.data.frame(results(dds, name = "time_24_vs_0", test = "Wald"))
time_24 <- tibble::rownames_to_column(time_24, "Geneid")
time_24$trt <- "24_100"
time_48<- as.data.frame(results(dds, name = "time_48_vs_0", test = "Wald"))
time_48 <- tibble::rownames_to_column(time_48, "Geneid")
time_48$trt <- "48_100"
time_72<- as.data.frame(results(dds, name = "time_72_vs_0", test = "Wald"))
time_72 <- tibble::rownames_to_column(time_72, "Geneid")
time_72$trt <- "72_100"
time_168<- as.data.frame(results(dds, name = "time_168_vs_0", test = "Wald"))
time_168 <- tibble::rownames_to_column(time_168, "Geneid")
time_168$trt <- "168_100"

time_0_50 <- as.data.frame(results(dds, name = "map_50_vs_100", test = "Wald"))
time_0_50 <- tibble::rownames_to_column(time_0_50, "Geneid")
time_0_50$trt <- "0_50"


time_3_map50 <- as.data.frame(results(dds, contrast=list(c("time_3_vs_0","map50.time3")), test = "Wald"))
time_3_map50 <- tibble::rownames_to_column(time_3_map50, "Geneid")
time_3_map50$trt <- "3_50"
time_24_map50 <- as.data.frame(results(dds, contrast=list(c("time_24_vs_0","map50.time24")), test = "Wald"))
time_24_map50 <- tibble::rownames_to_column(time_24_map50, "Geneid")
time_24_map50$trt <- "24_50"
time_48_map50 <- as.data.frame(results(dds, contrast=list(c("time_48_vs_0","map50.time48")), test = "Wald"))
time_48_map50 <- tibble::rownames_to_column(time_48_map50, "Geneid")
time_48_map50$trt <- "48_50"
time_72_map50 <- as.data.frame(results(dds, contrast=list(c("time_72_vs_0","map50.time72")), test = "Wald"))
time_72_map50 <- tibble::rownames_to_column(time_72_map50, "Geneid")
time_72_map50$trt <- "72_50"
time_168_map50 <- as.data.frame(results(dds, contrast=list(c("time_168_vs_0","map50.time168")), test = "Wald"))
time_168_map50 <- tibble::rownames_to_column(time_168_map50, "Geneid")
time_168_map50$trt <- "168_50"

de_time <- bind_rows(time_3, time_3_map50, time_24, time_24_map50,
                     time_48, time_48_map50, time_72, time_72_map50,
                     time_168, time_168_map50, time_0_50)

de_time_w_annotations <- left_join(de_time, annotations)

de_time %>%
  filter(trt == "3_100")%>%
  ggplot(., aes(baseMean, log2FoldChange))+
  #scale_x_continuous(trans='log10')+
  geom_point()

de_time %>%
  filter(trt == "3_100")%>%
  ggplot(., aes(baseMean))+
  geom_histogram()



#Get normalized cts. Then check some of the zeros to make sure they stayed zero
moo <- assay(dds)
vsd <- vst(dds, blind=TRUE)
moo2 <- assay(vsd)

#rld <- rlog(dds, blind=TRUE)
#moo3 <- assay(rld)

ntd <- normTransform(dds)
moo4 <- assay(ntd)

dds <- estimateSizeFactors(dds)
normalized_counts <- as.data.frame(counts(dds, normalized = T))

normalized_counts <- normalized_counts %>%
  mutate(Geneid = rownames(normalized_counts))



write.csv(de_time, file = "data/4w_de_cts_FC.csv")

write.csv(normalized_counts, file = "data/Normalized_cts_FC.csv")

de_time <- left_join(de_time, annotations)

#vmags <- as.list(read_csv("~/Documents/Projects/Wet-up/vMAG_no_duplicats.txt"))

temp <- de_time %>%
  filter(str_detect(kegg_id, "K03320")) %>%
  mutate(time = as.numeric(str_match(trt, "([0-9]+)_([0-9]+)")[,2]))

temp %>%
  ggplot(., aes(time, log2FoldChange, group = time))+
  geom_boxplot()




