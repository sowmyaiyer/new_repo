#awk '{ if ($6 == "-") {tss=$3} else if ($6 == "+") {tss=$2} printf("%s\t%d\t%d\t%s\t%s\t%s\n", $1,tss,tss,$4,$5,$6) }' ../txt/hg19_one_tr_per_gene.bed | bedtools slop -i stdin -g ~/gnearline/common_scripts/hg19.genome -r 150 -l 0 -s > ../txt/hg19_one_tr_per_gene.150downstream.bed
#awk '{ if ($6 == "-") { tss=$3} else if ($6 == "+") { tss=$2} printf("%s\t%d\t%d\t%s_%s_%d_%d\t%s\t%s\t%s\n",$1,tss,tss,$4,$1,$2,$3,$5,$6,$7)}' ../txt/hg19_snRNA_and_snoRNA.tmp.bed | bedtools slop -i stdin -g ~/gnearline/common_scripts/hg19.genome -r 150 -l 0 -s > ../txt/hg19_snRNA_and_snoRNA.tmp.150downstream.bed 
#cat ../txt/hg19_one_tr_per_gene.150downstream.bed ../txt/hg19_snRNA_and_snoRNA.tmp.150downstream.bed | sort -k4,4 > ../txt/hg19_one_tr_per_gene_with_snRNA_and_snoRNA.downstream150.bed 
#join -v 1 ../txt/hg19_one_tr_per_gene_with_snRNA_and_snoRNA.downstream150.bed ../txt/snRNAs_to_be_deleted.sorted.txt -1 4 -2 1 | sed 's/ /\t/g' | awk '{ print $2"\t"$3"\t"$4"\t"$1"\t"$5"\t"$6}' > ../txt/hg19_one_tr_per_gene_with_snRNA_and_snoRNA.downstream150.paralogs_removed.bed
for time in {"A","B","C","D","E","F","G","H"}
do
        echo """
#        bedtools intersect -nonamecheck -a ../txt/hg19_one_tr_per_gene_with_snRNA_and_snoRNA.downstream150.paralogs_removed.bed -b ../bowtie_out/112213flucap_${time}.sorted.bamtobed.bed -c -F 1 -s > ../txt/112213flucap_${time}.hg19_one_tr_per_gene_with_snRNA_and_snoRNA.downstream150.paralogs_removed.bed
#	awk '{ print \$4\"\\t\"\$7 }' ../txt/112213flucap_${time}.hg19_one_tr_per_gene_with_snRNA_and_snoRNA.downstream150.paralogs_removed.bed > ../txt/112213flucap_${time}.hg19_one_tr_per_gene_with_snRNA_and_snoRNA.downstream150.paralogs_removed.cap_counts_only.txt
#	join ../txt/112213flucap_${time}.hg19_one_tr_per_gene_with_snRNA_and_snoRNA.downstream150.paralogs_removed.cap_counts_only.txt ../txt/gene_mappings_unique_reads.${time}.txt  | sed 's/ /\t/g' > ../txt/capAndSnatchReadCounts.${time}.txt
#	echo "id multiplicity genenames" | sed 's/ /\t/g' > ../txt/112213flucap_${time}.multi_mapping_read_counts.proper.txt
#	sort -k3,3 ../txt/112213flucap_${time}.multi_mapping_read_counts.txt | bedtools groupby -i stdin -g 3 -c 1 -o count | awk '{ print \"line_\"NR\"\\t\"\$2\"\\t\"\$1}' | sort -k2,2nr >>  ../txt/112213flucap_${time}.multi_mapping_read_counts.proper.txt
	Rscript fluSnatch.R /home/si14w/gnearline/flu/txt/capAndSnatchReadCounts.${time}.txt  /home/si14w/gnearline/flu/txt/112213flucap_${time}.multi_mapping_read_counts.proper.txt ${time}
        """ > ../bsubFiles/capseq_readcounts_flu_snatch_contribution_${time}.hg19_one_tr_per_gene_with_snRNA_and_snoRNA.downstream150.paralogs_removed.bsub
done
