# This script measures the length between the three catalytic residues (H, C, D)
# The ULP protein structures were predicted using AlphaFold2 and the pdb files are available at https://doi.org/10.5061/dryad.1ns1rn97f
# talbe_of_cat-sites.csv is obtained from get_catalytic-site_position_from_MSA.py or see Table S4 

import os
import csv
from Bio import PDB
import pandas as pd

folder_path = './folder-name'
table_file_path = 'talbe_of_cat-sites.csv'
output_csv_path = 'output_name.csv'

position_dict = {}

with open(table_file_path, 'r') as csv_file:
    csv_reader = csv.reader(csv_file)
    # Skip the header row
    next(csv_reader)
    for row in csv_reader:
        name = row[0]
        h_position = int(row[1])
        d_position = int(row[2])
        c_position = int(row[3])
        position_dict[name] = (h_position, d_position, c_position)

result_dict = {}

pdb_parser = PDB.PDBParser(QUIET=True)
for name in position_dict:
    for pdb_filename in os.listdir(folder_path):
        if pdb_filename.endswith('.pdb'):
            pdb_path = os.path.join(folder_path, pdb_filename)

            # Parse the PDB file
            structure = pdb_parser.get_structure('protein', pdb_path)

            if pdb_filename.split('.pdb')[0] == name:
                model = structure[0]
                chain = model['A']
                residue1 = chain[position_dict[name][0]]
                residue2 = chain[position_dict[name][1]]
                residue3 = chain[position_dict[name][2]]

                atom1 = residue1['CA']
                atom2 = residue2['CA']
                atom3 = residue3['CA']

                distance12 = atom1 - atom2
                distance13 = atom1 - atom3
                distance23 = atom2 - atom3

                result_dict[name] = (distance12, distance13, distance23)



# saving the dataframe
df = pd.DataFrame.from_dict(result_dict, orient='index', columns=['H_D_Distance (Å)', 'H_C_Distance (Å)', 'D_C_Distance (Å)'])
df.to_csv(output_csv_path)
