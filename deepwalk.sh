#!/bin/bash

# FILE PARA DESCARGAR LOS PACKETES QUE SE USAN
python packages.py


f='deepwalk-embeddings/w10-80walks-40l-d128-'
b='.txt'
g='data/edge_lists_ints/'
X=("rattus_norvegicus" "inter_H-Y" "elegans" "coli" "drosop" "arab" "human" "mus_musculus" "yeast")

dir="deepwalk-embeddings"
if [[ ! -e $dir ]]; then
    mkdir $dir
fi

for w in {0..9}
do
	d=$f${X[$w]}$b
	deepwalk --format "edgelist" --input $g${X[$w]}$b --max-memory-data-size 0 --number-walks 80 --representation-size 128 --walk-length 40 --window-size 10 --seed 0 --workers 4 --output $d

done

python analyze_network.py
python process.py

