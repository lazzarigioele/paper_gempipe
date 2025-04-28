

# mamba install -c bioconda -c conda-forge 'carveme==1.6.2' 'python=3.9'
# python ~/cplex/python/setup.py install


cp -r ../gempipe/working/proteomes .
cp ../../reference/iYL1228_edited.xml .
cp ../../reference/medium_carveme.tsv .
cp ../../../helper_functions.py .


mkdir -p output
mkdir -p biolog_sims


for f in proteomes/*.faa
do
    echo "Processing $f";
    b=$(basename $f .faa);

    # run carveme:
    carve \
        -o output/$b.xml \
        -u gramneg \
        --fbc2 \
        --reference iYL1228_edited.xml \
        -g medium_carveme \
        -i medium_carveme \
        --mediadb medium_carveme.tsv \
        proteomes/$b.faa ;
    echo "Done processing $f $b";
done


for f in output/*.xml
do
    python << EOF
from helper_functions import *
simulate_biolog('$f', 'biolog_sims', seed=False, dataset='01_klebsiella')
EOF
done


