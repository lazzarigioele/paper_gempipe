
mkdir -p genomes_all

ncbi-genome-download \
     --no-cache \
     --metadata-table raw_ncbi_1598.txt \
     --retries 100 --parallel 10 \
     --output-folder genomes_all \
     --species-taxids 1598 \
     --formats assembly-stats,fasta \
     --section genbank \
     bacteria

mv raw_ncbi_1598.txt genomes_all/

