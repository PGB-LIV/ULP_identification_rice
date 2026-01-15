#!/bin/bash

# dowload dockers files from https://hub.docker.com/u/pegi3s
# run this script in a folder containing fasta files of ULP from yeast, Arabidopsis and RPRP rice

for x in *_Sc_At.fasta; 
do
A=${x%.*} # remove .fasta surfice; if only once ##should check the file name before
echo $A
docker run -it --rm  -v "path_to_folder":/data -w /data pegi3s/clustalomega -i $x -o "Sc_At"_"$A"_ClustalOmega.fasta -v;
docker run -it --rm  -v "path_to_folder":/data -w /data pegi3s/mafft bash -c "mafft --auto $x > "Sc_At"_"$A"_MAFFT.fasta"
docker run --rm  -v "path_to_folder":/data -w /data pegi3s/muscle -in $x -out "Sc_At"_"$A"_MUSCLE.fasta;
done
