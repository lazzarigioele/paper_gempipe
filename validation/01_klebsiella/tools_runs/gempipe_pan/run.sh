

cp -r ../../../01_klebsiella/genomes_renamed .
cp ../../reference_pan/KpSC_pan_metabolic_model_v2.xml .
cp ../../reference_pan/KpSC_pan_metabolic_model_v2_prots.faa .
cp ../../reference/medium_gempipe.json .


gempipe autopilot \
    -c 37 \
    -s neg \
    -g genomes_renamed/ \
    -b enterobacterales_odb10 \
    -o output/ \
    --verbose \
    --dbs ../../../../dbs/ \
    --buscoM 100% --ncontigs 10000 --N50 0 \
    -rm KpSC_pan_metabolic_model_v2.xml \
    -rp KpSC_pan_metabolic_model_v2_prots.faa \
    --media medium_gempipe.json \
    -rs KPN_SPONT \
    --minpanflux 0.9 \
    --minflux 0.6 \
    --biolog 
    
    
    
