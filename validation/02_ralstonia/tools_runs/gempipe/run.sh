

cp -r ../../genomes_renamed .
cp ../../reference/iRP1476_edited.xml .
cp ../../reference/iRP1476.faa .
cp ../../reference/medium_gempipe.json .
cp ../../../helper_functions.py .


mkdir -p biolog_sims


gempipe autopilot \
    -c 11 \
    -s neg \
    -g genomes_renamed/ \
    -b burkholderiales_odb10 \
    -o output/ \
    --verbose \
    --dbs ../../../../dbs/ \
    --buscoM 100% --ncontigs 10000 --N50 0 \
    -rm iRP1476_edited.xml \
    -rp iRP1476.faa \
    --media medium_gempipe.json \
    -rs spontaneous \
    --minpanflux 0.9 \
    --minflux 0.6 \
    --biolog \
    --mancor mancor.txt


# update simulations for pan-model:
rm output/biolog_panmodel.csv
python << EOF
from helper_functions import *
simulate_biolog('output/draft_panmodel.json', 'output', seed=False, dataset='02_ralstonia')
EOF
mv output/draft_panmodel.csv output/biolog_panmodel.csv


# update simulations for strain-specific models:
for f in output/strain_models_gf/*.json
do
    python << EOF
from helper_functions import *
simulate_biolog('$f', 'biolog_sims', seed=False, dataset='02_ralstonia')
EOF
done