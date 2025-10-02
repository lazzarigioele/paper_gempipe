

cp -r ../../proteomes_renamed .
cp ../../reference/medium_gempipe_rf.json .


gempipe autopilot \
    -c 9 \
    -s neg \
    -p proteomes_renamed/ \
    -b pseudomonas_odb12 \
    -o output/ \
    --verbose \
    --dbs ../../../../dbs/ \
    --buscoM 100% --ncontigs 10000 --N50 0 \
    --media medium_gempipe_rf.json \
    --minpanflux 0.9 \
    --minflux 0.6
    #--biolog 

