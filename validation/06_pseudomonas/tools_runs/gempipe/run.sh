

cp -r ../../proteomes_renamed .
cp ../../reference/iSD1509_M9_edited.xml .
cp ../../reference/iSD1509.faa .
cp ../../reference/medium_gempipe.json .


gempipe autopilot \
    -c 9 \
    -s neg \
    -p proteomes_renamed/ \
    -b pseudomonas_odb12 \
    -o output/ \
    --verbose \
    --dbs ../../../../dbs/ \
    --buscoM 100% --ncontigs 10000 --N50 0 \
    -rm iSD1509_M9_edited.xml \
    -rp iSD1509.faa \
    --media medium_gempipe.json \
    -rs PA_s0001 \
    --minpanflux 0.8 \
    --minflux 0.6
    #--biolog 
    #-rs 
    #--mancor mancor.txt 