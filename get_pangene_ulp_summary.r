# This script produces the Figure 1
# pangene_matrix_genes.tr.tab can be downloaded from https://zenodo.org/records/14772953 (Os4530.POR.tar.gz)

# load library
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggh4x)

# Pangene_version May 2024
pangene_table <- read_tsv("path_to_file/Os4530.POR.1/pangene_matrix_genes.tr.tab", col_names = T)
pangene_table_edited <- as.data.frame(apply(pangene_table,2, str_remove_all, "gene:")) %>% select(-18)
colnames(pangene_table_edited) <- c("Pan_id", "Nipponbare_merged", "OsAzu","OsCMeo","OsPr106","OsKeNa","OsARC", "OsZS97","OsN22","OsMH63","OsNaBo","OsLiXu","OsGoSa","OsLaMu","OsKYG", "OsIR64","OsLima")

Ulp_magic16_all_name <- read.csv("gene_with_ulp_catalyticSite_fullL_nonRedun_mod.csv")
pangene_ulp <- pangene_table_edited %>% filter_all( any_vars(. %in% Ulp_magic16_all_name$Gene))

#longer table with all nested values separated to multiple rows
only_panID_fullTable_sepLong <- pangene_ulp %>%
  separate_longer_delim(., 2:4, delim = ",") %>%
  separate_longer_delim(., 5, delim = ",") %>%
  separate_longer_delim(., 6, delim = ",") %>%
  separate_longer_delim(., 7:11, delim = ",") %>%
  separate_longer_delim(., 12:17, delim = ",")

# count the number of present of genes in each gene cluster
only_panID_fullTable_sepLong <- only_panID_fullTable_sepLong %>%  mutate(count_rice = rowSums(.[2:17] != '-'))

# make it long format
only_panID_fullTable_sepLong_if <- only_panID_fullTable_sepLong %>% gather(key = "Accession", value = "gene", -Pan_id, -count_rice)

# is ULP?
only_panID_fullTable_sepLong_if_classified <- only_panID_fullTable_sepLong_if %>%
  mutate(Is_ULP = if_else(gene %in% Ulp_magic16_all_name$Gene, "Yes", if_else(gene == "-","NA","No")))

# reorder
only_panID_fullTable_sepLong_if_classified$Is_ULP <- factor(only_panID_fullTable_sepLong_if_classified$Is_ULP, levels = c("Yes", "No", "NA"))

# rename Nipponbare_merged -> Osativa
only_panID_fullTable_sepLong_if_classified$Accession <- str_replace(only_panID_fullTable_sepLong_if_classified$Accession, "Nipponbare_merged","OsNip")

# Add grouping
rice_group <- read.csv("rice_grouping.csv")
only_panID_fullTable_sepLong_if_classified <- left_join(only_panID_fullTable_sepLong_if_classified, rice_group)

# reorder the grouping
only_panID_fullTable_sepLong_if_classified$Group <- factor(only_panID_fullTable_sepLong_if_classified$Group, levels = c("Tempjap rice","Tropjap rice","Aro rice","Aus rice","Indica rice"))


# Plot a summary of ULP across the rice genomes
color5 <- c("#ff4500","#ffa500","#ff6961","#800080","#ff69b4")
strip <- strip_themed(background_x = elem_list_rect(fill = color5))

p_1_8 <- ggplot(data = only_panID_fullTable_sepLong_if_classified %>% filter(count_rice <9), aes(x = Accession, y = reorder(Pan_id,desc(count_rice)), fill = Is_ULP)) +
  geom_tile(color = "black") +
  scale_fill_manual(values = c("red", "hotpink", "white"))+
  ylab("Pangene IDs")+
  xlab("Rice") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1),
        axis.text.y = element_text(size = 6),
        strip.text.x = element_text(color = "white", size = 0)) +
  facet_grid2(rows = vars(count_rice), cols = vars(Group), scales = "free", space = "free", strip = strip)

p_9_16 <- ggplot(data = only_panID_fullTable_sepLong_if_classified %>% filter(count_rice >8), aes(x = Accession, y = reorder(Pan_id,desc(count_rice)), fill = Is_ULP)) +
  geom_tile(color = "black") +
  scale_fill_manual(values = c("red", "hotpink","white"))+
  ylab("Pangene IDs")+
  xlab("Rice") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1),
        strip.text.x = element_text(color = "white", size = 0)) +
  facet_grid2(rows = vars(count_rice), cols = vars(Group), scales = "free", space = "free", strip = strip)

png(filename = "summary_ulp_pangene_ALL_ULP_corrected_PanGene_ver_May2024_reordered.png", height = 8.5, width = 9, units = "in", res = 300)
ggarrange(p_1_8, p_9_16, ncol = 2, common.legend = T)
dev.off()

