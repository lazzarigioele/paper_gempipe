

cp -r ../../proteomes_renamed .
cp ../../reference/Hirose2024_edited.xml .
cp ../../reference/iYH543.faa .
cp ../../reference/medium_gempipe.json .


gempipe autopilot \
    -c 17 \
    -s pos \
    -p proteomes_renamed/ \
    -b streptococcaceae_odb12 \
    -o output/ \
    --verbose \
    --dbs ../../../../dbs/ \
    --buscoM 100% --ncontigs 10000 --N50 0 \
    -rm Hirose2024_edited.xml \
    -rp iYH543.faa \
    --media medium_gempipe.json \
    --minpanflux 0.7 \
    --minflux 0.6
    #--biolog 
    #-rs 
    #--mancor mancor.txt 