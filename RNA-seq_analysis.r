# Download expression values (FPKM) by submitting the ULP gene names and search in https://plantrnadb.com/ricerna/ 
# FPKM.csv is a table of FPKM values with gene ID as columns and RNA-seq samples as rows. The last 8 columns are RNAseq data details: Tissue, Cultivar, Genotype, Treatment, Project, TotalReads, UniqueMappedRatio, ReleaseDate
# pan_id_MSU_model.csv is a table of 4 columns: Pan_id, MSU, count_rice, and count_domain
# note: count_rice = number of rice (out of 16) that has a gene in this pan-id
#       count_domain = number of rice (out of 16) that has a ULP gene in this pan-id
# example:
#|         Pan_id	       |     MSU	    |count_rice|count_domain|
#|Os4530.POR.1.pan0000812|LOC_Os01g02270|	8	       |    3       |


# load library
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggh4x)
library(patchwork)
library(gridExtra)
library(viridis)

panid <- read.csv("pan_id_MSU_model.csv") %>% separate_rows(MSU, sep = ",")
FPKM <- read.csv("FPKM.csv")
FPKM_value <- FPKM%>% select(Sample, grep("LOC", colnames(.)))

FPKM_value_selected_rice <- FPKM %>% filter(rowSums(.[,3:44] != 0) > 0) %>% 
  filter(str_detect(Cultivar, "Nipponbare|aponica|IR64|Minghui|Zhenshan|Azucena|N22|PR106|nivara|rufipogon")) %>% 
  filter(!str_detect(Cultivar, "x|X|Pandan|Tainung67|taipei|Minghui 86")) %>% 
  mutate(rice_group = case_when(str_detect(Cultivar, "nivara|rufipogon") ~ "1",
                                str_detect(Cultivar, "Azucena|Azuenca") ~ "3",
                                str_detect(Cultivar, "N22") ~ "4",
                                str_detect(Cultivar, "IR64|Minghui|Zhenshan|Azucena|PR106") ~ "5",
                                str_detect(Cultivar, "Nipponbare|aponica") ~ "2",
                                TRUE ~ "No")) 

FPKM_wild_rice <- FPKM %>% filter(str_detect(Cultivar, "nivara|rufipogon"))


all_FPKM_lf <- FPKM_value_selected_rice %>% select(-UniqueMappedRatio, -ReleaseDate, -TotalReads) %>% gather(., key = "gene_id", value = "FPKM", -Sample, -SampleName, -Cultivar, -Genotype, -Treatment, -Tissue, -rice_group, -Project) %>% 
  mutate(gene_id = gsub("S","s", .$gene_id)) %>% mutate(gene_id = gsub("G","g", .$gene_id)) %>% 
  left_join(., , by = c("gene_id"="MSU"))

all_FPKM_lf$rice_group <- factor(all_FPKM_lf$rice_group, levels = c("1","2","3","4","5"))



#unique(japonica$Treatment)
japonica_root<- all_FPKM_lf %>% filter(rice_group == "2", str_detect(Tissue, "root"))
japonica_root_uniq <- japonica_root %>% distinct(Sample) %>% mutate(order_plot = seq(1, nrow(.)))
japonica_root_mod <- japonica_root %>% left_join(japonica_root_uniq,by = "Sample")


p_japonica_root_1 <- ggplot(data = japonica_root_mod %>% filter(order_plot < 294), aes(x = order_plot, y = reorder(Pan_id,desc(count_domain)), fill = FPKM)) +
  geom_tile(color = "black") +
  ylab("Pangene IDs")+
  xlab("Sample") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text = element_text(colour = "black"),
        axis.text.y = element_text(size = 10),
        strip.text.x = element_text(color = "white", size = 0)) +
  scale_x_continuous(breaks=seq(0, 293, 10), expand = c(0, 0))+
  scale_fill_gradientn(colors = c("white","#E0B0FF","purple"), breaks  = c(0,5,10), limits = c(0,10),oob = scales::squish)+
  facet_grid2(rows = vars(count_domain), scales = "free", space = "free")

p_japonica_root_2 <- ggplot(data = japonica_root_mod%>% filter(order_plot > 293), aes(x = order_plot, y = reorder(Pan_id,desc(count_domain)), fill = FPKM)) +
  geom_tile(color = "black") +
  ylab("Pangene IDs")+
  xlab("Sample") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text = element_text(colour = "black"),
        axis.text.y = element_text(size = 10),
        strip.text.x = element_text(color = "white", size = 0)) +
  scale_x_continuous(breaks=seq(300, 586, 10), expand = c(0, 0))+
   scale_fill_gradientn(colors = c("white","#E0B0FF","purple"), breaks  = c(0,5,10), limits = c(0,10),oob = scales::squish)+
  facet_grid2(rows = vars(count_domain), scales = "free", space = "free")


png(filename = "FPKM_newsetULP_japonica_root.png", height = 15, width = 20, units = "in", res = 300)
p_japonica_root_1/p_japonica_root_2 + plot_layout(guides = "collect")
dev.off()

japonica_shoot<- all_FPKM_lf %>% filter(rice_group == "2", str_detect(Tissue, "shoot"))
japonica_shoot_uniq <- japonica_shoot %>% distinct(Sample) %>% mutate(order_plot = seq(1, nrow(.)))
japonica_shoot_mod <- japonica_shoot %>% left_join(japonica_shoot_uniq,by = "Sample")

p_japonica_shoot <- ggplot(data = japonica_shoot_mod , aes(x = order_plot, y = reorder(Pan_id,desc(count_domain)), fill = FPKM)) +
  geom_tile(color = "black") +
  ylab("Pangene IDs")+
  xlab("Sample") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text = element_text(colour = "black"),
        axis.text.y = element_text(size = 10),
        strip.text.x = element_text(color = "white", size = 0)) +
  scale_x_continuous(breaks=seq(0, 380, 10), expand = c(0, 0))+
  scale_fill_gradientn(colors = c("white","#E0B0FF","purple"), breaks  = c(0,5,10), limits = c(0,10),oob = scales::squish)+
  facet_grid2(rows = vars(count_domain), scales = "free", space = "free")


png(filename = "FPKM_newsetULP_japonica_shoot.png", height = 7.5, width = 20, units = "in", res = 300)
p_japonica_shoot
dev.off()

japonica_leaf<- all_FPKM_lf %>% filter(rice_group == "2", str_detect(Tissue, "leaf"))
japonica_leaf_uniq <- japonica_leaf %>% distinct(Sample) %>% mutate(order_plot = seq(1, nrow(.)))
japonica_leaf_mod <- japonica_leaf %>% left_join(japonica_leaf_uniq,by = "Sample")


wild_rice <- all_FPKM_lf %>% filter(rice_group == "1")
wild_rice_uniq <- wild_rice %>% distinct(Sample) %>% mutate(order_plot = seq(1, nrow(.)))
wild_rice_mod <- wild_rice %>% left_join(wild_rice_uniq,by = "Sample")


p_wild_rice <- ggplot(data = wild_rice_mod , aes(x = order_plot, y = reorder(Pan_id,desc(count_domain)), fill = FPKM)) +
  geom_tile(color = "black") +
  ylab("Pangene IDs")+
  xlab("Sample") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text = element_text(colour = "black"),
        axis.text.y = element_text(size = 10),
        strip.text.x = element_text(color = "white", size = 0)) +
  scale_x_continuous(breaks=seq(0, 12, 1), expand = c(0, 0))+
  scale_fill_gradientn(colors = c("white","#E0B0FF","purple"), breaks  = c(0,5,10), limits = c(0,10),oob = scales::squish)+
  facet_grid2(rows = vars(count_domain), scales = "free", space = "free")


png(filename = "FPKM_newsetULP_wild_rice.png", height = 7.5, width = 10, units = "in", res = 300)
p_wild_rice
dev.off()

library(ComplexHeatmap)
library(circlize)

#2486 RNA-seq
selected_data <- rbind(wild_rice,japonica_shoot, japonica_leaf, japonica_root, magic_rice) 

FPKM_value_selected_forPlot <- FPKM_value_selected_rice %>% filter(Sample %in% selected_data$Sample) %>% select(1:45, 48) %>% left_join(., selected_data %>% select(Sample, rice_group) %>% distinct())


colnames(FPKM_value_selected_forPlot) <- gsub("G","g", colnames(FPKM_value_selected_forPlot))
colnames(FPKM_value_selected_forPlot) <- gsub("S","s", colnames(FPKM_value_selected_forPlot))


col_order <- data.frame(MSU = colnames(FPKM_value_selected_forPlot)[3:44]) %>% left_join(., )

#1759
FPKM_value_selected_forPlot_japonica <- FPKM_value_selected_forPlot %>% filter(rice_group == 2)
input <- FPKM_value_selected_forPlot_japonica[, -c(2, 45:47)] %>% column_to_rownames(var = "sample")
input_scale <- scale(input)

anno_tissue <- data.frame(Tissue = FPKM_value_selected_forPlot_japonica$Tissue) %>% 
                          mutate(reduce_tissue = case_when(str_detect(Tissue, "leaf") ~ "leaf",
                                                    str_detect(Tissue, "shoot") ~ "shoot",
                                                    str_detect(Tissue, "root") ~ "root",
                                                    TRUE ~ "no"))

row_ha <- rowAnnotation(Tissue = anno_tissue$reduce_tissue, col = list(Tissue = c("leaf" = "darkgreen", "shoot"="brown", "root"="orange")))

png("Heatmap_Zcolumn_japonica.png", width = 9, height = 11, units = "in", res = 300)
p_jap <- Heatmap(input_scale,
        name = "z score\nby column",
        column_split =  col_order$count_domain, 
        cluster_columns = T, 
        cluster_rows = T, 
        show_row_names = F,
        cluster_row_slices = F, 
        cluster_column_slices = F,
        column_gap = unit(4, "mm"),
        col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        border = T,
        left_annotation = row_ha)
p_jap
dev.off()


ht <- draw(p_jap)

order_tissue <- rownames(input_scale)[row_order(ht)]

order_tissue_info <-FPKM_value_selected_forPlot %>% select(1:2,45:46) %>% column_to_rownames(var = "sample")
order_tissue_info <- order_tissue_info[order_tissue,]

write.csv(order_tissue_info, "order_tissue_info.csv")


#631
FPKM_value_selected_forPlot_japonica_onlyControl <- FPKM_value_selected_forPlot %>% filter(rice_group == 2, str_detect(Treatment, "control|Control|mock|--"))
input <- FPKM_value_selected_forPlot_japonica_onlyControl[, -c(2, 45:47)] %>% column_to_rownames(var = "sample")
input_scale <- scale(input)

anno_tissue <- data.frame(Sample = FPKM_value_selected_forPlot_japonica_onlyControl$sample, Tissue = FPKM_value_selected_forPlot_japonica_onlyControl$Tissue) %>% 
                          mutate(reduce_tissue = case_when(str_detect(Tissue, "shoot|whole|stem") ~ "shoot",
                                                           str_detect(Tissue, "leaf") ~ "leaf",
                                                           str_detect(Tissue, "root") ~ "root",
                                                    TRUE ~ "no"))


row_ha <- rowAnnotation(Tissue = anno_tissue$reduce_tissue, col = list(Tissue = c("leaf" = "darkgreen", "shoot"="brown", "root"="orange")))

png("Heatmap_Zcolumn_japonica_onlyControl.png", width = 9, height = 11, units = "in", res = 300)
p_jap_con <- Heatmap(input_scale,
        name = "z score\nby column",
        column_split =  col_order$count_domain, 
        cluster_columns = T, 
        cluster_rows = T, 
        show_row_names = F,
        cluster_row_slices = F, 
        cluster_column_slices = F,
        column_gap = unit(4, "mm"),
        col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        border = T,
        left_annotation = row_ha)
p_jap_con
dev.off()

ht <- draw(p_jap_con)

order_tissue <- rownames(input_scale)[row_order(ht)]

order_tissue_info <-FPKM_value_selected_forPlot %>% select(1:2,45:46) %>% column_to_rownames(var = "sample")
order_tissue_info <- order_tissue_info[order_tissue,]

write.csv(order_tissue_info, "order_tissue_info_onlyControl.csv")

# summary plot
anno_tissue2 <- data.frame(Sample = FPKM_value_selected_forPlot_japonica_onlyControl$sample, Tissue = FPKM_value_selected_forPlot_japonica_onlyControl$Tissue) %>% 
                          mutate(reduce_tissue = case_when(str_detect(Tissue, "shoot|whole|stem") ~ "shoot",
                                                           str_detect(Tissue, "leaf") ~ "leaf",
                                                           str_detect(Tissue, "root") ~ "root",
                                                    TRUE ~ "no"))
#use mean of the FPKM
sum_plot <- as.data.frame(input_scale) %>% rownames_to_column(var = "Sample") %>% left_join(., anno_tissue2) %>% select(-Tissue) %>% group_by(reduce_tissue)%>% summarise(across(-Sample, mean, na.rm = TRUE)) %>% column_to_rownames(var = "reduce_tissue")

sum_plot_rearranged <- sum_plot %>% select(c("LOC_Os03g24990","LOC_Os08g33280","LOC_Os10g24954","LOC_Os11g12780","LOC_Os02g16240","LOC_Os12g24880","LOC_Os01g02270","LOC_Os06g03180","LOC_Os05g10270","LOC_Os11g12500","LOC_Os03g41780","LOC_Os05g34520","LOC_Os05g37550", "LOC_Os07g32090","LOC_Os04g30860","LOC_Os06g13750","LOC_Os12g01290","LOC_Os07g12990", "LOC_Os11g42610","LOC_Os04g54670","LOC_Os07g18300","LOC_Os03g24960","LOC_Os03g24980","LOC_Os01g33520","LOC_Os03g25000","LOC_Os04g49394","LOC_Os04g28590","LOC_Os01g25370","LOC_Os06g29310","LOC_Os12g36580","LOC_Os03g29630","LOC_Os11g01180","LOC_Os10g06910","LOC_Os09g08440","LOC_Os02g27280","LOC_Os10g33450","LOC_Os09g12480","LOC_Os11g10780","LOC_Os12g41380","LOC_Os03g22400","LOC_Os05g11770","LOC_Os06g28030"))


png("Heatmap_Zcolumn_japonica_onlyControl_sumbyTissue.png", width = 9, height = 3, units = "in", res = 300)
p_sum <- Heatmap(sum_plot_rearranged,
        name = "mean\nz score",
        column_split =  c(rep(1,5), rep(3,3), rep(5,6), rep(6,2),8, rep(9,3),11 ,rep(12,4),13,14, rep(15,3), rep(16,12) ), 
        cluster_columns = F, 
        cluster_rows = F, 
        show_row_names = T,
        row_names_side = "left",
        cluster_row_slices = F, 
        cluster_column_slices = F,
        column_gap = unit(4, "mm"),
        col=colorRamp2(c(-0.5, 0,0.5), c("blue", "white", "red")),
        border = T)
p_sum
dev.off()


# to scale by RNA-seq of the same gene (by row)
input <- FPKM_value_selected_forPlot_japonica_onlyControl[, -c(2, 45:47)] %>% column_to_rownames(var = "sample")
input_scale_row <- scale(t(input))

row_ha <- rowAnnotation(Tissue = anno_tissue2$reduce_tissue, col = list(Tissue = c("leaf" = "darkgreen", "shoot"="brown", "root"="orange")))

png("Heatmap_Zrow_japonica_onlyControl.png", width = 9, height = 11, units = "in", res = 300)
p_jap_con_row <- Heatmap(t(input_scale_row),
        name = "z score\nby row",
        column_split =  col_order$count_domain, 
        cluster_columns = T, 
        cluster_rows = T, 
        show_row_names = F,
        cluster_row_slices = F, 
        cluster_column_slices = F,
        column_gap = unit(4, "mm"),
        col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        border = T,
        left_annotation = row_ha)
p_jap_con_row

japonica <- all_FPKM_lf %>% filter(rice_group == "2")
japonica$gene_id <- factor(japonica$gene_id, levels = c("LOC_Os08g33280", "LOC_Os03g24990", "LOC_Os11g12780", "LOC_Os02g16240", "LOC_Os10g24954","LOC_Os01g02270", "LOC_Os06g03180","LOC_Os12g24880","LOC_Os07g32090","LOC_Os11g12500", "LOC_Os05g10270", "LOC_Os03g41780", "LOC_Os05g37550", "LOC_Os05g34520","LOC_Os04g30860","LOC_Os06g13750","LOC_Os12g01290","LOC_Os11g42610","LOC_Os07g12990", "LOC_Os04g54670","LOC_Os07g18300","LOC_Os03g24980","LOC_Os03g24960","LOC_Os01g33520","LOC_Os03g25000","LOC_Os04g49394","LOC_Os04g28590","LOC_Os06g29310","LOC_Os01g25370","LOC_Os12g36580","LOC_Os11g10780", "LOC_Os03g22400", "LOC_Os12g41380", "LOC_Os05g11770","LOC_Os09g12480", "LOC_Os10g33450", "LOC_Os09g08440", "LOC_Os10g06910", "LOC_Os02g27280","LOC_Os06g28030", "LOC_Os11g01180","LOC_Os03g29630"))



library(ggtext) 
japonica_2 <- all_FPKM_lf %>% filter(rice_group == "2")
japonica_2$gene_id <- factor(japonica_2$gene_id, levels = c("LOC_Os03g24990","LOC_Os08g33280","LOC_Os10g24954","LOC_Os11g12780","LOC_Os02g16240","LOC_Os12g24880","LOC_Os01g02270","LOC_Os06g03180","LOC_Os05g10270","LOC_Os11g12500","LOC_Os03g41780","LOC_Os05g34520","LOC_Os05g37550", "LOC_Os07g32090","LOC_Os04g30860","LOC_Os06g13750","LOC_Os12g01290","LOC_Os07g12990", "LOC_Os11g42610","LOC_Os04g54670","LOC_Os07g18300","LOC_Os03g24960","LOC_Os03g24980","LOC_Os01g33520","LOC_Os03g25000","LOC_Os04g49394","LOC_Os04g28590","LOC_Os01g25370","LOC_Os06g29310","LOC_Os12g36580","LOC_Os03g29630","LOC_Os11g01180","LOC_Os10g06910","LOC_Os09g08440","LOC_Os02g27280","LOC_Os10g33450","LOC_Os09g12480","LOC_Os11g10780","LOC_Os12g41380","LOC_Os03g22400","LOC_Os05g11770","LOC_Os06g28030"))


# Salt
japonica_salt <- japonica_2 %>% filter(str_detect(Treatment, "Salt|salt|NaCl"))
japonica_salt_control <- japonica_2 %>% filter(Project %in% japonica_salt$Project) %>%
  filter(!str_detect(Treatment, "cold|Drought|Heat|Cold|Highlight|drought|--")) %>% 
  mutate(Treatment_short = case_when(str_detect(Treatment, "control|Control") ~ "control",
                                     str_detect(Treatment, "Salt|salt|NaCl") ~"salt",
                                     TRUE ~ "no"))
p_japonica_salt_control<-ggplot(data = japonica_salt_control, aes(x=gene_id, y = FPKM, fill = Treatment_short)) +
  geom_point(alpha = 0.5, size = 0.4, position = "jitter")+
  geom_boxplot(width=0.8, outliers = F) +
  facet_grid2(cols = vars(count_domain),  scales = "free", space = "free")+
  theme_bw() +
  scale_fill_manual(values = c("#4e91fd", "#FA5F55")) +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        panel.spacing.x = unit(0.5, "lines"),
        strip.text.x = element_text(size = 14)) +
  ggtitle("Control vs Salt stress") +
  stat_compare_means( label = "p.signif", label.y = 130, hide.ns = T, colour = "red", na.rm = T)

#Drought
japonica_drought <- japonica_2 %>% filter(str_detect(Treatment, "Drought|drought|dehy"))
japonica_drought_control <- japonica_2 %>% filter(Project %in% japonica_drought$Project) %>%
  filter(!str_detect(Treatment, "cold|cadmium|Heat|Cold|Highlight|--|Flood|Salt|salt|NaCl|ABA|JA|Day"),
         !str_detect(SampleName, "cold|cadmium|Heat|Cold|Highlight|--|Flood|Salt|salt|NaCl|ABA|JA|Cadmium")) %>% 
  mutate(Treatment_short = case_when(str_detect(Treatment, "control|Control") ~ "control",
                                     str_detect(Treatment, "Drought|drought|Dehydra|Osmotic") ~"drought",
                                     TRUE ~ "no"))
p_japonica_drought_control<-ggplot(data = japonica_drought_control, aes(x=gene_id, y = FPKM, fill = Treatment_short)) +
  geom_point(alpha = 0.5, size = 0.4, position = "jitter")+
  geom_boxplot(width=0.8, outliers = F) +
  facet_grid2(cols = vars(count_domain),  scales = "free", space = "free")+
  theme_bw() +
  scale_fill_manual(values = c("#4e91fd", "#FA5F55")) +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        panel.spacing.x = unit(0.5, "lines"),
        strip.text.x = element_text(size = 14))  +
  ggtitle("Control vs Drought stress")+
  stat_compare_means( label = "p.signif", label.y = 80, hide.ns = T, colour = "red", na.rm = T)

# Heat
japonica_heat <- japonica_2 %>% filter(str_detect(Treatment, "heat|Heat"))
japonica_heat_control <- japonica_2 %>% filter(Project %in% japonica_heat$Project) %>%
  filter(!str_detect(Treatment, "cold|Drought|Cold|Highlight|drought|--|Salt|salt")) %>% 
  mutate(Treatment_short = case_when(str_detect(Treatment, "control|Control") ~ "control",
                                     str_detect(Treatment, "heat|Heat") ~"heat",
                                     TRUE ~ "no"))
p_japonica_heat_control<-ggplot(data = japonica_heat_control, aes(x=gene_id, y = FPKM, fill = Treatment_short)) +
  geom_point(alpha = 0.5, size = 0.4, position = "jitter")+
  geom_boxplot(width=0.8, outliers = F) +
  facet_grid2(cols = vars(count_domain),  scales = "free", space = "free")+
  theme_bw() +
  scale_fill_manual(values = c("#4e91fd", "#FA5F55")) +
  theme(axis.text.x = element_blank(),
        axis.title  = element_blank(),
        panel.spacing.x = unit(0.5, "lines"),
        strip.text.x = element_text(size = 14))  +
  ggtitle("Control vs Heat stress") +
  stat_compare_means( label = "p.signif", label.y = 70, hide.ns = T, colour = "red", na.rm = T)


# ABA 
japonica_ABA <- japonica_2 %>% filter(str_detect(Treatment, "ABA"))
japonica_ABA_control <- japonica_2 %>% filter(Project %in% japonica_ABA$Project)%>%
  filter(!str_detect(Treatment, "cold|Drought|Cold|Osmotic|drought|--|Flood|cadmium|ethylene|JA|GA")) %>% 
  filter(!str_detect(SampleName, "cold|Drought|Cold|Osmotic|drought|Development|Flood|cadmium|ethylene|JA|GA|Cadmium")) %>% 
  mutate(Treatment_short = case_when(str_detect(Treatment, "mock|control") ~ "control",
                                     str_detect(Treatment, "ABA") ~"ABA",
                                     TRUE ~ "no"))

japonica_ABA_control$Treatment_short <- factor(japonica_ABA_control$Treatment_short, levels = c("control", "ABA"))

p_japonica_ABA_control<-ggplot(data = japonica_ABA_control, aes(x=gene_id, y = FPKM, fill = Treatment_short)) +
  geom_point(alpha = 0.5, size = 0.4, position = "jitter")+
  geom_boxplot(width=0.8, outliers = F) +
  facet_grid2(cols = vars(count_domain),  scales = "free", space = "free")+
  theme_bw() +
  scale_fill_manual(values = c("#4e91fd", "#FA5F55")) +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        panel.spacing.x = unit(0.5, "lines"),
        strip.text.x = element_text(size = 14)) +
  ggtitle("Control vs ABA treatment") +
  stat_compare_means( label = "p.signif", label.y = 78, hide.ns = T, colour = "red", na.rm = T)


# Infection
japonica_infection <- japonica_2 %>% filter(str_detect(Treatment, "infection|infected"))
japonica_infection_control <- japonica_2 %>% filter(Project %in% japonica_infection$Project) %>% 
  mutate(Treatment_short = case_when(str_detect(Treatment, "mock|control") ~ "control",
                                     str_detect(Treatment, "infection|infected|transcected|treatted|Meloidogyne") ~"infection",
                                     TRUE ~ "no"))
p_japonica_infection_control<-ggplot(data = japonica_infection_control, aes(x=gene_id, y = FPKM, fill = Treatment_short)) +
  geom_point(alpha = 0.5, size = 0.4, position = "jitter")+
  geom_boxplot(width=0.8, outliers = F) +
  facet_grid2(cols = vars(count_domain),  scales = "free", space = "free")+
  theme_bw() +
  scale_fill_manual(values = c("#4e91fd", "#FA5F55")) +
  theme(axis.text.x = element_text(angle = 90, colour = "black", size = 16, hjust = 0.5, vjust = 0.5),
        axis.title = element_blank(),
        plot.title = element_markdown(),
        panel.spacing.x = unit(0.5, "lines"),
        strip.text.x = element_text(size = 14)) +
  labs(title = "Control vs nematode/fungi/virus infection") +
  stat_compare_means( label = "p.signif", label.y = 80, hide.ns = T, colour = "red", na.rm = T)

png("4stress_FPKM_ULP_japonica_newOrder_all.png", width = 14, height = 15, units = "in", res = 400)
grid.arrange(patchworkGrob(p_japonica_salt_control / p_japonica_drought_control/p_japonica_heat_control/p_japonica_ABA_control/p_japonica_infection_control), left = "FPKM") 
dev.off()
