

mkdir -p gempipe
cd gempipe


gempipe recon \
    -c 36 \
    -s pos \
    -g ../genomes_tf/ \
    -b lactobacillales_odb10 \
    -o output_recon/ \
    --verbose \
    --dbs ../../dbs/ \
    -rm ../reference/Lreuteri_530_fixed.json \
    -rp ../reference/proteome_supplemented.faa \
    --sbml \
    --nofig \
    --buscoM 3% --buscoF 2% --ncontigs 240 --N50 19000 
