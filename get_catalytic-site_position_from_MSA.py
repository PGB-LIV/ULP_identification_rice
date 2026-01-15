# This script identify the position of the 3 catalystic residues (H, C, D) on ULP proteins

import os
from Bio import AlignIO
from Bio import SeqIO
import pandas as pd

# path to a folder containg MSA (ClustalOmega, MAFFT, MUSCLE) files in fasta format
folder_path = "D:/Os4530.POR.1/Os4530.POR.1/ULP_clusters/"

# Can also use the length to compare
# ScULP1 - last char is the catalytic site 
# H
seq_1 = "MSVEVDKHRNTLQYHKKNPYSPLFSPISTYRCYPRVLNNPSESRRSASFSGIYKKRTNTSRFNYLNDRRVLSMEESMKDGSDRASKAGFIGGIRETLWNSGKYLWHTFVKNEPRNFDGSEVEASGNSDVESRSSGSRSSDVPYGLRENYSSDTRKHKFDTSTWALPNKRRRIESEGVGTPSTSPISSLASQKSNCDSDNSITFSRDPFGWNKWKTSAIGSNSENNTSDQKNSYDRRQYGTAFIRKKKVAKQNINNTKLVSRAQSEEVTYLRQIFNGEYKVPKILKEERERQLKLMDMDKEKDTGLKKSIIDLTEKIKTILIENNKNRLQTRNENDDDLVFVKEKKISSLERKHKDYLNQKLKFDRSILEFEKDFKRYNEILNERKKIQEDLKKKKEQLAKKKLVPELNEKDDDQVQKALASRENTQLMNRDNIEITVRDFKTLAPRRWLNDTIIEFFMKYIEKSTPNTVAFNSFFYTNLSERGYQGVRRWMKRKKTQIDKLDKIFTPINLNQSH" 
# D
seq_2 = seq_1 + "WALGIIDLKKKTIGYVD"
# C
seq_3 = seq_2 + "SLSNGPNAMSFAILTDLQKYVMEESKHTIGEDFDLIHLDCPQQPNGYDC"


for filename in os.listdir(folder_path):
  if filename.endswith(("MUSCLE.fasta", "ClustalOmega.fasta", "MAFFT.fasta")):
      
    file_path = os.path.join(folder_path, filename)
    file_name = filename.split(".")[0]
    print("Started: "+ file_name)
    
    # read the alignment
    alignment = AlignIO.read(file_path, "fasta")
    number_of_seqs = len(alignment)
    id_name = []
    
    
    # variables to store values
    pos_H = -1
    pos_D = -1
    pos_C = -1
    
    # list with gene names
    for seq_nucleotides in alignment:
        gene_name = seq_nucleotides.id
        id_name.append(gene_name.split(' ')[0])
    
    # Search for the position of the catalytic sites in the alignment
    for i in range(number_of_seqs):
      if id_name[i] == "ScULP1": # use as reference

        #this could collapse into one loop
        for random_1 in range(100, 2500):
          seq_random_1 = str(alignment[i].seq[0:random_1]).replace("-", "")
          if seq_random_1 == seq_1:  
            #print("found H at position", random_1) 
            pos_H = random_1
            break
                
        for random_2 in range(100, 2500):
          seq_random_2 = str(alignment[i].seq[0:random_2]).replace("-", "")
          if seq_random_2 == seq_2:  
            #print("found D at position", random_2) 
            pos_D = random_2
            break
            
        for random_3 in range(100, 2500):
          seq_random_3 = str(alignment[i].seq[0:random_3]).replace("-", "")
          if seq_random_3 == seq_3:  
            #print("found C at position", random_3)
            pos_C = random_3
            break

            
    #H
    pos_interest_1 = pos_H - 1 #minus one
    length_to_pos_1 = []
    chars_at_pos_1 = []
    #D
    pos_interest_2 = pos_D - 1 #minus one
    length_to_pos_2 = []
    chars_at_pos_2 = []
    #C
    pos_interest_3 = pos_C - 1 #minus one
    length_to_pos_3 = []
    chars_at_pos_3 = []
    
    
    for i in range(0,number_of_seqs):
        seq_H = alignment[i].seq[0: pos_interest_1].replace("-", "")
        length_to_pos_1.append(len(seq_H)+1) #add 1, seq starts from 1
        chars_at_pos_1.append(alignment[i].seq[pos_interest_1])
    
        seq_D = alignment[i].seq[0: pos_interest_2].replace("-", "")
        length_to_pos_2.append(len(seq_D) + 1)  # add 1, seq starts from 1
        chars_at_pos_2.append(alignment[i].seq[pos_interest_2])
    
        seq_C = alignment[i].seq[0: pos_interest_3].replace("-", "")
        length_to_pos_3.append(len(seq_C) + 1)  # add 1, seq starts from 1
        chars_at_pos_3.append(alignment[i].seq[pos_interest_3])
    
    # dictionary of lists
    dict = {'Gene': id_name,
            'AA1': chars_at_pos_1, 'Pos1': length_to_pos_1,
            'AA2': chars_at_pos_2, 'Pos2': length_to_pos_2,
            'AA3': chars_at_pos_3, 'Pos3': length_to_pos_3}
    
    df = pd.DataFrame(dict)
    
    # saving the dataframe
    df.to_csv("./catalytic_position_results/" + file_name +'_catalytic-domain_positions.csv')
    
    print("Finished: " + file_name)
