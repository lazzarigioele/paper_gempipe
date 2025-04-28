

cp -r ../../genomes_renamed .
cp ../../reference/iYL1228_edited.xml .
cp ../../reference/iYL1228.faa .
cp ../../reference/medium_gempipe.json .


gempipe autopilot \
    -c 24 \
    -s neg \
    -g genomes_renamed/ \
    -b enterobacterales_odb10 \
    -o output/ \
    --verbose \
    --dbs ../../../../dbs/ \
    --buscoM 100% --ncontigs 10000 --N50 0 \
    -rm iYL1228_edited.xml \
    -rp iYL1228.faa \
    --media medium_gempipe.json \
    -rs KPN_SPONT \
    --minpanflux 0.9 \
    --minflux 0.6 \
    --biolog \
    --mancor mancor.txt 
