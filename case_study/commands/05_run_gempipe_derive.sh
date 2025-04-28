

mkdir -p gempipe
cd gempipe


gempipe derive --verbose \
    -c 36 \
    -o output_derive/ \
    -im draft_panmodel_edit.json \
    -ip output_recon/pam.csv \
    -ir report_updated.csv \
    -m cdm_reuteri.json \
    --minflux 0.5 \
    --aux \
    --cnps \
    --biosynth 0.8