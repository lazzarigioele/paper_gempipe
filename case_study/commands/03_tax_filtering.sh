

mkdir -p fastani
mkdir -p fastani/input_fna

cp genomes_all/genbank/bacteria/*/*.fna fastani/input_fna  

ANIclustermap \
	-i fastani/input_fna \
	-o fastani \
	--mode fastani \
	--thread_num 36 \
	--overwrite \
	--fig_width 10 \
	--fig_height 10 \
	--annotation 

