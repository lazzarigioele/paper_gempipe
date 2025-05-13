

mkdir -p fastani


ANIclustermap \
	-i genomes_nomags \
	-o fastani \
	--mode fastani \
	--thread_num 36 \
	--overwrite \
	--fig_width 10 \
	--fig_height 10 \
	--annotation 

