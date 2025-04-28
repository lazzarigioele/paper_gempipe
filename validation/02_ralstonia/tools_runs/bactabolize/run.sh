

# mamba install -c kelwyres -c bioconda -c conda-forge 'bactabolize==1.0.3' 'python=3.9'
# python ~/cplex/python/setup.py install 
# pip install gempipe
# mamba install depinfo==1.7.0


cp -r ../gempipe/working/proteomes .
cp ../../reference/iRP1476_edited.xml .
cp ../../reference/iRP1476.faa .
cp ../../reference/iRP1476.fna .
cp ../../reference/medium_bactabolize.json .
cp ../../../helper_functions.py .


mkdir -p output
mkdir -p gapfilled
mkdir -p biolog_sims
mkdir -p genbanks


#for f in proteomes/*.faa
#do
#    python << EOF
#from helper_functions import *
#convert_fasta_to_genbank_with_cds('$f', 'genbanks')
#EOF
#done


## copy and rename medium definition:
#cp medium_bactabolize.json /opt/Progs/miniconda3/envs/bactabolize/lib/python3.9/site-#packages/bactabolize/data/media_definitions/medium_bactabolize_02.json


#for f in genbanks/*.gb
#do
#    b=$(basename $f);
#    bactabolize draft_model \
#        --assembly_fp $f \
#        --no_reannotation \
#        --ref_proteins_fp iRP1476.faa \
#        --ref_genes_fp iRP1476.fna \
#        --ref_model_fp iRP1476_edited.xml \
#        --biomass_reaction_id BIOMASS \
#        --media_type medium_bactabolize_02 \
#        --atmosphere_type aerobic \
#        --output_fp $b ;
#    # --min_coverage 25 and --min_pident 80 are default
#    b=$(basename $f .gb);
#    mv $b* output/ ;
#done


for f in output/*.json
do {
    python << EOF
from helper_functions import *
gapfilling_for_bactabolize('iRP1476_edited.xml', '$f', 'gapfilled', dataset='02_ralstonia')
EOF
} & done; wait


for f in gapfilled/*.json
do {
    python << EOF
from helper_functions import *
simulate_biolog('$f', 'biolog_sims', seed=False, dataset='02_ralstonia', tool='bactabolize')
EOF
} & done; wait