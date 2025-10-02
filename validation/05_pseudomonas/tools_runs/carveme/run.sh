

# mamba install -c bioconda -c conda-forge 'carveme==1.6.2' 'python=3.9'
# python ~/cplex/python/setup.py install


cp -r ../../proteomes_renamed .
cp ../../reference/iSD1509_edited.xml .
cp ../../reference/medium_carveme.tsv .


mkdir -p output


for f in proteomes_renamed/*.faa
do {
    echo "Processing $f";
    b=$(basename $f .faa);

    # run carveme:
    carve \
        -o output/$b.xml \
        -u gramneg \
        --fbc2 \
        --reference iSD1509_edited.xml \
        -g medium_carveme \
        -i medium_carveme \
        --mediadb medium_carveme.tsv \
        proteomes_renamed/$b.faa ;
    echo "Done processing $f $b";
} & done; wait


