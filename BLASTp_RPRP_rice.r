# Download protein annotation files (fasta) of RPRP rice from https://ftp.gramene.org/oryza/PanOryza/fasta/ 
# Install BLAST programs locally: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
# The input file is the reference ULP protein sequences from yeast and Arabidopsis in fasta format

# Create the database by using the fasta files 
# note that the file names are not updated, see the latest names at https://ftp.gramene.org/oryza/PanOryza/fasta/ and https://panoryza.org/node/6
# Run this once for each rice, using command line (bash)
#makeblastdb -in oryza_sativair64.evd.protein.fasta -dbtype prot -out oryza_sativair64/oryza_sativair64
#makeblastdb -in oryza_sativaazucena.evd.protein.fasta -dbtype prot -out oryza_sativaazucena/oryza_sativaazucena
#makeblastdb -in oryza_sativa132424.evd.protein.fasta -dbtype prot -out oryza_sativa132424/oryza_sativa132424
#makeblastdb -in oryza_sativa132278.evd.protein.fasta -dbtype prot -out oryza_sativa132278/oryza_sativa132278
#makeblastdb -in oryza_sativa128077.evd.protein.fasta -dbtype prot -out oryza_sativa128077/oryza_sativa128077
#makeblastdb -in oryza_sativa127742.evd.protein.fasta -dbtype prot -out oryza_sativa127742/oryza_sativa127742
#makeblastdb -in oryza_sativa127652.evd.protein.fasta -dbtype prot -out oryza_sativa127652/oryza_sativa127652
#makeblastdb -in oryza_sativa127564.evd.protein.fasta -dbtype prot -out oryza_sativa127564/oryza_sativa127564
#makeblastdb -in oryza_sativa127518.evd.protein.fasta -dbtype prot -out oryza_sativa127518/oryza_sativa127518
#makeblastdb -in oryza_sativa125827.evd.protein.fasta -dbtype prot -out oryza_sativa125827/oryza_sativa125827
#makeblastdb -in oryza_sativa125619.evd.protein.fasta -dbtype prot -out oryza_sativa125619/oryza_sativa125619
#makeblastdb -in oryza_sativa117425.evd.protein.fasta -dbtype prot -out oryza_sativa117425/oryza_sativa117425
#makeblastdb -in oryza_aus.evd.protein.fasta -dbtype prot -out oryza_aus/oryza_aus
#makeblastdb -in oryza_mh63.evd.protein.fasta -dbtype prot -out oryza_mh63/oryza_mh63
#makeblastdb -in oryza_zs97.evd.protein.fasta -dbtype prot -out oryza_zs97/oryza_zs97
#makeblastdb -in oryza_sativa.evd.protein.fasta -dbtype prot -out oryza_sativa/oryza_sativa


library(tidyverse)
library(ggplot2)

genome_database <- c("oryza_sativair64","oryza_sativaazucena","oryza_sativa132424","oryza_sativa132278","oryza_sativa128077","oryza_sativa127742","oryza_sativa127652","oryza_sativa127564","oryza_sativa127518","oryza_sativa125827","oryza_sativa125619","oryza_sativa117425","oryza_sativa","oryza_aus", "oryza_mh63", "oryza_zs97")

#"oryza_aus", "oryza_mh63", "oryza_zs97" = diferent gene name system; do it separately

y_lab_list <- list()

for (i in 1:12) {
  
  blast_db = paste("path_to_file/proteins/",genome_database[i],"/",genome_database[i], sep = "")
  input = "./input_file.fasta"
  #evalue = 1e-6
  format = 6
  colnames <- c("qseqid",
                "sseqid",
                "pident",
                "length",
                "mismatch",
                "gapopen",
                "qstart",
                "qend",
                "sstart",
                "send",
                "e_value",
                "bitscore")
  
  blast_out <- system2(command = "blastp", 
                       args = c("-db", blast_db, 
                                "-query", input, 
                                "-outfmt", format),
                       wait = TRUE,
                       stdout = TRUE) %>%
    as_tibble() %>% 
    separate(col = value, 
             into = colnames,
             sep = "\t",
             convert = TRUE)
  
  blast_out_filtered <- blast_out %>% 
    separate(., sseqid, into = c("geneID","geneModel"), sep = "([.?:])") %>%
    filter(geneModel == "01") %>% select(qseqid, geneID, pident, length,e_value) %>%
    group_by(qseqid) 
  
  # modify output table: get the chromosome
  blast_out_filtered <-blast_out_filtered %>% mutate(chr1 = case_when(str_detect(qseqid,"Os01g|Os01t")~"1", 
                                                                      str_detect(qseqid,"Os02g|Os02t")~ "2",
                                                                      str_detect(qseqid,"Os03g|Os03t")~"3", 
                                                                      str_detect(qseqid,"Os04g|Os04t")~ "4",
                                                                      str_detect(qseqid,"Os05g|Os05t")~"5", 
                                                                      str_detect(qseqid,"Os06g|Os06t")~ "6",
                                                                      str_detect(qseqid,"Os07g|Os07t")~"7", 
                                                                      str_detect(qseqid,"Os08g|Os08t")~ "8",
                                                                      str_detect(qseqid,"Os09g|Os09t")~"9", 
                                                                      str_detect(qseqid,"Os10g|Os10t")~ "10",
                                                                      str_detect(qseqid,"Os11g|Os11t")~"11", 
                                                                      str_detect(qseqid,"Os12g|Os12t")~ "12",
                                                                      TRUE ~ "no"),
                                                     chr2 = case_when(str_detect(geneID,"_01g")~"1", 
                                                                      str_detect(geneID,"_02g")~ "2",
                                                                      str_detect(geneID,"_03g")~"3", 
                                                                      str_detect(geneID,"_04g")~ "4",
                                                                      str_detect(geneID,"_05g")~"5", 
                                                                      str_detect(geneID,"_06g")~ "6",
                                                                      str_detect(geneID,"_07g")~"7", 
                                                                      str_detect(geneID,"_08g")~ "8",
                                                                      str_detect(geneID,"_09g")~"9", 
                                                                      str_detect(geneID,"_10g")~ "10",
                                                                      str_detect(geneID,"_11g")~"11", 
                                                                      str_detect(geneID,"_12g")~ "12",
                                                                      TRUE ~ "no")) %>% 
    mutate(same_chr = if_else(chr1==chr2, "Yes","No"))
  
  # create and save figures
  y_lab <- blast_out_filtered %>% group_by(qseqid) %>% 
    count(same_chr) %>% spread(key = same_chr, value = n) %>% 
    mutate(color = if_else(Yes >= 1, "#0E7A34", "black", "black"))
  
  y_lab_list[[i]] <- y_lab # store output in the list
  
  temp_plot <- ggplot(data = blast_out_filtered, aes(x = geneID, y = qseqid, color = e_value, size = pident,shape = same_chr, alpha=same_chr)) +
    geom_point()+ 
    theme_bw()+ theme(axis.title = element_blank(), 
                      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"), 
                      axis.text.y = element_text(color = y_lab$color)) +
    scale_color_gradientn(colors = c("red", "blue"), limits = c(0,10^-6)) +
    scale_alpha_discrete(range = c(0.5, 1))+
    ggtitle(paste("BLASTP (top 3 by e_value): ", genome_database[i]), 
            "Green genes and triangel: hit to the same chromosome, Red: significant (e_value < 1e-6)")+
    guides(pident = guide_legend(order = 1), 
           e_value = guide_legend(order = 2),
           same_chr = guide_legend(order = 3))
  
  ggsave(temp_plot, file=paste("Blast_output/",genome_database[i],".png",sep = ""), height = 8, width = 17, units = "in", dpi = 300)
  
  
  y_lab <- blast_out_filtered %>% group_by(qseqid) %>% count(same_chr) %>% spread(key = same_chr, value = n) %>% mutate( sum = if_else(Yes >= 1, 1, 0, 0))
  colnames(y_lab) <- c("qseqid","No","Yes",genome_database[i] )
  
  #export table
  write.csv(blast_out_filtered, paste("Blast_output/blastP_out_",genome_database[i],".csv",sep = ""))
  write.csv(y_lab, paste("Blast_output/chr_location/chr_location_",genome_database[i],".csv",sep = ""))
}

#"oryza_aus", "oryza_mh63", "oryza_zs97" = diferent gene name system; do it separately
# different separation of gene name;Osaus.11G006060_01  see below;


for (i in 13:16) {
blast_db = paste("path_to_file/proteins/",genome_database[i],"/",genome_database[i], sep = "")

input = "./Mix-input_file.fasta"
#evalue = 1e-6
format = 6
colnames <- c("qseqid",
              "sseqid",
              "pident",
              "length",
              "mismatch",
              "gapopen",
              "qstart",
              "qend",
              "sstart",
              "send",
              "e_value",
              "bitscore")

blast_out <- system2(command = "blastp", 
                     args = c("-db", blast_db, 
                              "-query", input, 
                              "-outfmt", format),
                     wait = TRUE,
                     stdout = TRUE) %>%
  as_tibble() %>% 
  separate(col = value, 
           into = colnames,
           sep = "\t",
           convert = TRUE)

# select top 3 (least e-values)
blast_out_filtered <- blast_out %>% 
  separate(., sseqid, into = c("geneID","geneModel"), sep = "_") %>%
  filter(geneModel == "01") %>% select(qseqid, geneID, pident, length,e_value) %>%
  group_by(qseqid) %>% slice_min(.,n = 3, order_by = e_value)

# modify output table: get the chromosome
blast_out_filtered <-blast_out_filtered %>% mutate(chr1 = case_when(str_detect(qseqid,"Os01g|Os01t")~"1", 
                                                                    str_detect(qseqid,"Os02g|Os02t")~ "2",
                                                                    str_detect(qseqid,"Os03g|Os03t")~"3", 
                                                                    str_detect(qseqid,"Os04g|Os04t")~ "4",
                                                                    str_detect(qseqid,"Os05g|Os05t")~"5", 
                                                                    str_detect(qseqid,"Os06g|Os06t")~ "6",
                                                                    str_detect(qseqid,"Os07g|Os07t")~"7", 
                                                                    str_detect(qseqid,"Os08g|Os08t")~ "8",
                                                                    str_detect(qseqid,"Os09g|Os09t")~"9", 
                                                                    str_detect(qseqid,"Os10g|Os10t")~ "10",
                                                                    str_detect(qseqid,"Os11g|Os11t")~"11", 
                                                                    str_detect(qseqid,"Os12g|Os12t")~ "12",
                                                                    TRUE ~ "no"),
                                                   chr2 = case_when(str_detect(geneID,"01G")~"1", #notice capital G
                                                                    str_detect(geneID,"02G")~ "2",
                                                                    str_detect(geneID,"03G")~"3", 
                                                                    str_detect(geneID,"04G")~ "4",
                                                                    str_detect(geneID,"05G")~"5", 
                                                                    str_detect(geneID,"06G")~ "6",
                                                                    str_detect(geneID,"07G")~"7", 
                                                                    str_detect(geneID,"08G")~ "8",
                                                                    str_detect(geneID,"09G")~"9", 
                                                                    str_detect(geneID,"10G")~ "10",
                                                                    str_detect(geneID,"11G")~"11", 
                                                                    str_detect(geneID,"12G")~ "12",
                                                                    TRUE ~ "no")) %>% mutate(same_chr = if_else(chr1==chr2, "Yes","No"))

# create and save figures
y_lab <- blast_out_filtered %>% group_by(qseqid) %>% count(same_chr) %>% spread(key = same_chr, value = n) %>% mutate(color = if_else(Yes >= 1, "#0E7A34", "black", "black"))


temp_plot <- ggplot(data = blast_out_filtered, aes(x = geneID, y = qseqid, color = e_value, size = pident,shape = same_chr, alpha=same_chr)) +
  geom_point()+ 
  theme_bw()+ theme(axis.title = element_blank(), 
                    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"), 
                    axis.text.y = element_text(color = y_lab$color)) +
  scale_color_gradientn(colors = c("red", "blue"), limits = c(0,10^-6)) +
  scale_alpha_discrete(range = c(0.5, 1))+
  ggtitle(paste("BLASTP (top 3 by e_value): ", genome_database[i]), "Green genes and triangel: hit to the same chromosome, Red: significant (e_value < 1e-6)")+
  guides(pident = guide_legend(order = 1), 
         e_value = guide_legend(order = 2),
         same_chr = guide_legend(order = 3))

ggsave(temp_plot, file=paste("Blast_output/",genome_database[i],".png",sep = ""), height = 8, width = 17, units = "in", dpi = 300)


y_lab <- blast_out_filtered %>% group_by(qseqid) %>% count(same_chr) %>% spread(key = same_chr, value = n) %>% mutate( sum = if_else(Yes >= 1, 1, 0, 0))
colnames(y_lab) <- c("qseqid","No","Yes",genome_database[i] )

#export table
write.csv(blast_out_filtered, paste("Blast_output/blastP_out_",genome_database[i],".csv",sep = ""))
write.csv(y_lab, paste("Blast_output/chr_location/chr_location_",genome_database[i],".csv",sep = ""))
}

