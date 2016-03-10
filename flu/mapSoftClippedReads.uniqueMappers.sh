for time in {"A","B","C","D","E","F","G","H"}
do
	echo """
	samtools view ../bowtie_out/112213flucap_${time}.softclip.sorted.paralogs_removed.bam | awk '{ print \$1}' | sort | uniq -c | awk '{ if (\$1 == 1) print \$2 }' | sort > ../txt/112213flucap_${time}.unique_mapping_reads.txt
	samtools view ../bowtie_out/112213flucap_${time}.softclip.sorted.paralogs_removed.bam | sort -k1,1 > ../bowtie_out/112213flucap_${time}.softclip.sorted.paralogs_removed.sam
	join  ../bowtie_out/112213flucap_${time}.softclip.sorted.paralogs_removed.sam  ../txt/112213flucap_${time}.unique_mapping_reads.txt -1 1 -2 1 | sed 's/ /\t/g' | awk '{ print \$1\"\\t\"\$3}'| sort -k2,2 | bedtools groupby -g 2 -c 1 -o count > ../txt/gene_mappings_unique_reads.incomplete.${time}.txt
	join ../txt/hg19_one_tr_per_gene_and_snRNA_snoRNA.paralogs_removed.genenames.txt ../txt/gene_mappings_unique_reads.incomplete.${time}.txt -1 1 -2 1 -a 1 | awk '{ if (NF == 1 ) print \$1\"\\t\"0; else print \$1\"\\t\"\$2}' > ../txt/gene_mappings_unique_reads.${time}.txt
	samtools view ../bowtie_out/112213flucap_${time}.softclip.sorted.paralogs_removed.bam | awk '{ print \$1}' | sort | uniq -c | awk '{ if (\$1 > 1 && \$1 <= 20 ) print \$2 }' | sort > ../txt/112213flucap_${time}.multi_mapping_reads.txt
	join  ../bowtie_out/112213flucap_${time}.softclip.sorted.paralogs_removed.sam  ../txt/112213flucap_${time}.multi_mapping_reads.txt -1 1 -2 1 | sed 's/ /\\t/g' | sort -k1,1 | bedtools groupby -g 1 -c 3,3 -o count,distinct > ../txt/112213flucap_${time}.multi_mapping_read_counts.txt """ > ../bsubFiles/preProcess_em.${time}.bsub
done
