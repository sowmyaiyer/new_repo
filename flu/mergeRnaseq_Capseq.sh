for time in {"6h","12h","24h","48h","mock12h"}
do
	echo $time
	awk '{ print $4"\t"$8 }' ../txt/112213flucap_${time}.readCounts_150bp.rpm.bed | sort -k1,1 | awk '{ print $1}' > ../txt/delete.${time}.geneList.capseq.txt
	awk '{ print $4"\t"$8 }' ../txt/112213flucap_${time}.readCounts_150bp.rpm.bed | sort -k1,1  > ../txt/112213flucap_${time}.readCounts_150bp.rpm.sorted_by_gene.txt
	awk '{ if (NR > 1) print $1"\t"$6 }' ../txt/04SEP15.A549_B59_totalRNA_${time}.isoforms.results | sort -k1,1 | awk '{ print $1}' > ../txt/delete.${time}.geneList.rnaseq.txt
	awk '{ if (NR > 1) print $1"\t"$6 }' ../txt/04SEP15.A549_B59_totalRNA_${time}.isoforms.results | sort -k1,1 > ../txt/04SEP15.A549_B59_totalRNA_${time}.isoforms.sorted_by_gene.txt
	echo "diff  ../txt/delete.${time}.geneList.capseq.txt ../txt/delete.${time}.geneList.rnaseq.txt"
	diff ../txt/delete.${time}.geneList.capseq.txt ../txt/delete.${time}.geneList.rnaseq.txt
	join ../txt/112213flucap_${time}.readCounts_150bp.rpm.sorted_by_gene.txt ../txt/04SEP15.A549_B59_totalRNA_${time}.isoforms.sorted_by_gene.txt -1 1 -2 1 | sed 's/ /\t/g' > ../txt/capseq_vs_rnaseq.${time}.txt
	Rscript plotCapseqVsRnaseq.R $HOME/gnearline/flu/txt/capseq_vs_rnaseq.${time}.txt $HOME/gnearline/flu/txt/capseq_vs_rnaseq.${time}.percentiles.txt $HOME/gnearline/flu/txt/capseq_vs_rnaseq.${time}.pdf ${time}
done
#join ../txt/delete.${time}.geneList.rnaseq.txt ../txt/hg19_refseq-2.sorted.bed -1 1 -2 1 -a 1 > ucsc_genesymbol_table.txt
#paste ../txt/capseq_vs_rnaseq.6h.txt ../txt/capseq_vs_rnaseq.12h.txt ../txt/capseq_vs_rnaseq.24h.txt ../txt/capseq_vs_rnaseq.48h.txt ../txt/capseq_vs_rnaseq.mock12h.txt > ../txt/capseq_vs_rnaseq.all.txt
#awk '{ if (($3 > 10) && ($6 > 10) && ($9 > 10) && ($12 > 10) && ($15 > 10)) print $1"\t"$2"\t"$3"\t"}' ../txt/capseq_vs_rnaseq.all.txt > ../txt/capseq_vs_rnaseq.all.filtered_10tpm.txt
