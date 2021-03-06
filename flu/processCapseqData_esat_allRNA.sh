#awk '{ if ($6 == "-") {start=$3-150;end=$3} else if ($6 == "+") {start=$2;end=$2+150} printf("%s\t%d\t%d\t%s\t%s\t%s\n", $1,start,end,$4,$5,$6) }' ../txt/hg19_ucsc_genes_and_flu.bed > ../txt/hg19_ucsc_genes_and_flu.5prime_plus150.allRNA.bed
#sort -k4,4 hg19_snRNA.bed > hg19_snRNA.sorted_by_genename.bed
#sort -k1,1 hg19_ens_to_genename.bed > hg19_ens_to_genename.sorted_by_ens_genename.bed
#join hg19_snRNA.sorted_by_genename.bed hg19_ens_to_genename.sorted_by_ens_genename.bed -1 4 -2 1| awk -F" " '{ print $2"\t"$3"\t"$4"\t"$NF"\t"0"\t"$6}' | sort -k1,1 -k2,2n > ../txt/hg19_snRNA.sorted.bed
#awk '{ if ($6 == "-") {start=$3-150;end=$3} else if ($6 == "+") {start=$2;end=$2+150} printf("%s\t%d\t%d\t%s\t%s\t%s\n", $1,start,end,$4,$5,$6) }' ../txt/hg19_snRNA.sorted.bed >> ../txt/hg19_ucsc_genes_and_flu.5prime_plus150.allRNA.bed

#sort -k4,4 hg19_snoRNA.bed > hg19_snoRNA.sorted_by_genename.bed
#join hg19_snoRNA.sorted_by_genename.bed hg19_ens_to_genename.sorted_by_ens_genename.bed -1 4 -2 1| awk -F" " '{ print $2"\t"$3"\t"$4"\t"$NF"\t"0"\t"$6}' | sort -k1,1 -k2,2n > ../txt/hg19_snoRNA.sorted.bed
#awk '{ if ($6 == "-") {start=$3-150;end=$3} else if ($6 == "+") {start=$2;end=$2+150} printf("%s\t%d\t%d\t%s\t%s\t%s\n", $1,start,end,$4,$5,$6) }' ../txt/hg19_snoRNA.sorted.bed >> ../txt/hg19_ucsc_genes_and_flu.5prime_plus150.allRNA.bed

awk 'BEGIN{OFS="\t"}{ print $6,$11,$12}' /share/data/umw_biocore/genome_data/rsem_correction/refseq_hg19/truechr_hg19_ucsc_refseq_final.bed | head


for time in {"A","B","C","D","E","F","G","H"}
do
	echo """
	java -jar ./esat.v0.1_09.18.15_10.06.jar -in  ../bowtie_out/112213flucap_${time}.sorted.bam -annotations -geneMapping ../txt/knownIsoforms_hg19_and_flu.txt -out ../esat_out/esat_out_112213flucap_${time}.txt
	bedtools coverage -nonamecheck -a ../txt/hg19_ucsc_genes_and_flu.5prime_plus150.allRNA.bed -b ../bowtie_out/112213flucap_${time}.sorted.bamtobed.bed -counts -F 1 -s > ../txt/112213flucap_${time}.readCounts_150bp.allRNA.bed
	Rscript normalizeReadCounts.R $HOME/gnearline/flu/txt/112213flucap_${time}.readCounts_150bp.allRNA.bed $HOME/gnearline/flu/txt/112213flucap_${time}.readCounts_150bp.rpm.allRNA.bed
	""" > ../bsubFiles/capseq_readcounts_flu_${time}.bsub
done
