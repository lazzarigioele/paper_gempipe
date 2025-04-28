
# conda install gapseq==1.4 python=3.9 mamba
# python ~/cplex/python/setup.py install
# pip install gempipe


cp -r ../gempipe/working/proteomes .
cp ../../reference/medium_gapseq.csv .
cp ../../../helper_functions.py .


mkdir -p output
mkdir -p biolog_sims


for f in proteomes/*.faa
do {
    b=$(basename $f .faa);
    mkdir -p output/${b}/ ; 
    cp $f output/${b}/${b}.faa ;
    gapseq doall -K 1 output/${b}/${b}.faa medium_gapseq.csv ;  
    mv ${b}* output/${b}/ ;
    python << EOF
from helper_functions import *
simulate_biolog('output/${b}/${b}.xml', 'biolog_sims', seed=True, dataset='03_pseudomonas', tool='gapseq')
EOF
} & done; wait