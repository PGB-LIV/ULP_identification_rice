# Download protein annotation files (fasta) of wild rice from https://oryza-ensembl.gramene.org/species.html
# Install BLAST programs locally: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
# The input file is the reference ULP protein sequences from yeast and Arabidopsis in fasta format (https://doi.org/10.5061/dryad.1ns1rn97f)

# Create the database by using the fasta file
# Run this once for each rice, using command line (bash)
#makeblastdb -in Leersia_perrieri.Lperr_V1.4.pep.all.fa -dbtype prot -out Leersia_perrieri/Leersia_perrieri 
#makeblastdb -in Oryza_barthii.O.barthii_v1.pep.all.fa -dbtype prot -out Oryza_barthii/Oryza_barthii
#makeblastdb -in Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa -dbtype prot -out Oryza_brachyantha/Oryza_brachyantha
#makeblastdb -in Oryza_glaberrima.Oryza_glaberrima_V1.pep.all.fa -dbtype prot -out Oryza_glaberrima/Oryza_glaberrima
#makeblastdb -in Oryza_glumipatula.Oryza_glumaepatula_v1.5.pep.all.fa -dbtype prot -out Oryza_glumipatula/Oryza_glumipatula
#makeblastdb -in Oryza_longistaminata.O_longistaminata_v1.0.pep.all.fa -dbtype prot -out Oryza_longistaminata/Oryza_longistaminata
#makeblastdb -in Oryza_nivara.Oryza_nivara_v1.0.pep.all.fa -dbtype prot -out Oryza_nivara/Oryza_nivara
#makeblastdb -in Oryza_punctata.Oryza_punctata_v1.2.pep.all.fa -dbtype prot -out Oryza_punctata/Oryza_punctata
#makeblastdb -in Oryza_rufipogon.OR_W1943.pep.all.fa -dbtype prot -out Oryza_rufipogon/Oryza_rufipogon

# create BLast_output_all_hits and BLast_output_all_hits/chr_location  folders
system2(command = "mkdir", args = "BLast_output_all_hits" )
system2(command = "mkdir", args = "BLast_output_all_hits/chr_location" )

library(tidyverse)
library(ggplot2)

genome_database <- c("Leersia_perrieri","Oryza_barthii","Oryza_glaberrima",	"Oryza_glumipatula","Oryza_nivara","Oryza_punctata","Oryza_rufipogon", "Oryza_brachyantha",	"Oryza_longistaminata")

y_lab_list <- list()

for (i in 1:10) {
  
  blast_db = paste("./protein_databases/",genome_database[i],"/",genome_database[i], sep = "")
  input = "./input.fasta"
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
    filter(geneModel == "1") %>% select(qseqid, geneID, pident, length,e_value) %>%
    group_by(qseqid) 
  
  
  # modify output table: get the chromosome number 
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
                                                     chr2 = case_when(str_detect(geneID,"01G")~"1", 
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
                                                                      TRUE ~ "no")) %>% mutate(same_chr = if_else(chr1==chr2, "Yes","No"),
                                                      genome = paste(genome_database[i])
           )
  
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
    ggtitle(paste("BLASTP (all hits): ", genome_database[i]), 
            "Green genes and triangel: hit to the same chromosome, Red: significant (e_value < 1e-6)")+
    guides(pident = guide_legend(order = 1), 
           e_value = guide_legend(order = 2),
           same_chr = guide_legend(order = 3))
  
  ggsave(temp_plot, file=paste("BLast_output_all_hits/",genome_database[i],".png",sep = ""), height = 8, width = 25, units = "in", dpi = 300)
  
  
  y_lab <- blast_out_filtered %>% group_by(qseqid) %>% count(same_chr) %>% spread(key = same_chr, value = n) %>% mutate( sum = if_else(Yes >= 1, 1, 0, 0))
  colnames(y_lab) <- c("qseqid","No","Yes",genome_database[i] )
  
  #export table
  write.csv(blast_out_filtered, paste("BLast_output_all_hits/blastP_out_",genome_database[i],".csv",sep = ""))
  write.csv(y_lab, paste("BLast_output_all_hits/chr_location/chr_location_",genome_database[i],".csv",sep = ""))
}

