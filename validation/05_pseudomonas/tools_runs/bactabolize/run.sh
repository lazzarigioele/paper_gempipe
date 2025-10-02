

# mamba install -c kelwyres -c bioconda -c conda-forge 'bactabolize==1.0.3' 'python=3.9'
# python ~/cplex/python/setup.py install 
# pip install gempipe
# mamba install depinfo==1.7.0

 
cp -r ../../proteomes_renamed .
cp ../../reference/iSD1509_edited.xml .
cp ../../reference/iSD1509.faa .
cp ../../reference/iSD1509.fna .
cp ../../reference/medium_bactabolize.json .
cp ../../../helper_functions.py .


mkdir -p output
mkdir -p gapfilled
mkdir -p genbanks


for f in proteomes_renamed/*.faa
do { 
    python << EOF
from helper_functions import *
convert_fasta_to_genbank_with_cds('$f', 'genbanks')
EOF
} & done; wait


## copy and rename medium definition:
cp medium_bactabolize.json /opt/Progs/miniconda3/envs/bactabolize/lib/python3.9/site-packages/bactabolize/data/media_definitions/medium_bactabolize_05.json


for f in genbanks/*.gb
do {
    b=$(basename $f);
    bactabolize draft_model \
        --assembly_fp $f \
        --no_reannotation \
        --ref_proteins_fp iSD1509.faa \
        --ref_genes_fp iSD1509.fna \
        --ref_model_fp iSD1509_edited.xml \
        --biomass_reaction_id BIOMASS_PA14_v27M \
        --media_type medium_bactabolize_05 \
        --atmosphere_type aerobic \
        --output_fp $b 
    b=$(basename $f .gb);
    mv $b* output/ ;
} & done; wait


for f in output/*.json
do {
    python << EOF
from helper_functions import *
gapfilling_for_bactabolize('iSD1509_edited.xml', '$f', 'gapfilled', dataset='05_streptococcus')
EOF
} & done; wait

