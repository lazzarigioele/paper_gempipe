

cp -r ../../proteomes_renamed .
cp ../../reference/medium_gempipe_rf.json .


gempipe autopilot \
    -c 17 \
    -s pos \
    -p proteomes_renamed/ \
    -b streptococcaceae_odb12 \
    -o output/ \
    --verbose \
    --dbs ../../../../dbs/ \
    --buscoM 100% --ncontigs 10000 --N50 0 \
    --media medium_gempipe_rf.json \
    --minpanflux 0.7 \
    --minflux 0.6
    #--biolog 

