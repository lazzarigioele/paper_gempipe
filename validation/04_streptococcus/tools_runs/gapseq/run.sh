
# conda install gapseq==1.4 python=3.9 mamba
# python ~/cplex/python/setup.py install
# pip install gempipe


cp -r ../../proteomes_renamed .
cp ../../reference/medium_gapseq_dedup.csv .
cp ../../../helper_functions.py .


mkdir -p output


for f in proteomes_renamed/*.faa
do {
    b=$(basename $f .faa);
    mkdir -p output/${b}/ ; 
    cp $f output/${b}/${b}.faa ;
    gapseq doall -K 1 output/${b}/${b}.faa medium_gapseq_dedup.csv ;  
    mv ${b}* output/${b}/ ;
} & done; wait