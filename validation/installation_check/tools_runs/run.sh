

cp -r ../../01_klebsiella/genomes_renamed .
cp ../reference/KpSC_pan_metabolic_model_v2.json .
cp ../reference/KpSC_pan_metabolic_model_v2_prots.faa .


gempipe recon --verbose \
    -c 8 \
    -g genomes_renamed/ \
    -o output1/ 
    

gempipe recon --verbose \
    -c 8 \
    -g genomes_renamed/ \
    -o output2/ \
    -rm KpSC_pan_metabolic_model_v2.json \
    -rp KpSC_pan_metabolic_model_v2_prots.faa 
    
    
gempipe autopilot --verbose \
    -c 8 \
    -g genomes_renamed/ \
    -o output3/ \
    -rm KpSC_pan_metabolic_model_v2.json \
    -rp KpSC_pan_metabolic_model_v2_prots.faa 
