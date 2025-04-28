

cp -r ../../genomes_renamed .
cp ../../reference/iJN1463_edited.xml .
cp ../../reference/iJN1463.faa .
cp ../../reference/medium_gempipe.json .


gempipe autopilot \
    -c 24 \
    -s neg \
    -g genomes_renamed/ \
    -b pseudomonadales_odb10 \
    -o output/ \
    --verbose \
    --dbs ../../../../dbs/ \
    --buscoM 100% --ncontigs 10000 --N50 0 \
    -rm iJN1463_edited.xml \
    -rp iJN1463.faa \
    --media medium_gempipe.json \
    -rs PP_s0001 \
    --minpanflux 0.9 \
    --minflux 0.6 \
    --biolog \
    --mancor mancor.txt 