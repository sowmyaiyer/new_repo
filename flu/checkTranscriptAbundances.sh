for time in {"6h","12h","24h","48h","mock12h"}
do
	join ../txt/04SEP15.A549_B59_totalRNA_${time}.isoforms.sorted_by_gene.txt ../txt/hg19_refseq-2.sorted.bed | sed 's/ /\t/g' | sort -k3,3 | bedtools groupby -i stdin -g 3 -c 2 -o collapse | sed 's/,/\t/g' > ../txt/04SEP15.A549_B59_totalRNA_${time}.transcript_abundances.txt
done
