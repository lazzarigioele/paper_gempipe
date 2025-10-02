

mkdir -p cocoremover
mkdir -p cocoremover/input
mkdir -p cocoremover/output

for line in $(cat tables/list_MAGs.txt); do
    
    cp genomes_all/genbank/bacteria/${line}/${line}*.fna cocoremover/input/${line}.fna
    cocoremover -c 36 -i cocoremover/input/${line}.fna -t 1598 -d cocoremover/cocoremover.db -o cocoremover/output

done

