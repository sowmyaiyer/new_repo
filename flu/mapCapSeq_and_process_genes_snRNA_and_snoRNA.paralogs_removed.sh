for time in {"A","B","C","D","E","F","G","H"}
do
        echo """	
	module load bowtie2/2-2.1.0
	bowtie2 -x  ../bowtie_out/hg19_snRNA_and_snoRNA.sequences.150.paralogs_removed.fa -U /project/umw_garberlab/narayanan/112213flucap/rawData_toFluAndtoHuman/${time}_Nocodes.fastq.gz -S ../bowtie_out/112213flucap_${time}.hg19_one_tr_per_gene_and_snRNA_snoRNA.paralogs_removed.sam
	samtools view -bS -F 4 -o ../bowtie_out/112213flucap_${time}.hg19_one_tr_per_gene_and_snRNA_snoRNA.paralogs_removed.bam ../bowtie_out/112213flucap_${time}.hg19_one_tr_per_gene_and_snRNA_snoRNA.paralogs_removed.sam
	samtools sort -T ../bowtie_out/112213flucap_${time}.sorted -o ../bowtie_out/112213flucap_${time}.sorted.hg19_one_tr_per_gene_and_snRNA_snoRNA.paralogs_removed.bam ../bowtie_out/112213flucap_${time}.hg19_one_tr_per_gene_and_snRNA_snoRNA.paralogs_removed.bam
	samtools index ../bowtie_out/112213flucap_${time}.sorted.hg19_one_tr_per_gene_and_snRNA_snoRNA.paralogs_removed.bam
	bedtools bamtobed -i ../bowtie_out/112213flucap_${time}.sorted.hg19_one_tr_per_gene_and_snRNA_snoRNA.paralogs_removed.bam > ../bowtie_out/112213flucap_${time}.sorted.bamtobed.hg19_one_tr_per_gene_and_snRNA_snoRNA.paralogs_removed.bed
	bedtools groupby -i ../bowtie_out/112213flucap_${time}.sorted.bamtobed.hg19_one_tr_per_gene_and_snRNA_snoRNA.paralogs_removed.bed -g 1 -c 4 -o count > ../txt/capSeqReadCounts_paralogs_removed_tss_plus_150.${time}.txt
	rm  ../bowtie_out/112213flucap_${time}.hg19_one_tr_per_gene_and_snRNA_snoRNA.paralogs_removed.sam  ../bowtie_out/112213flucap_${time}.hg19_one_tr_per_gene_and_snRNA_snoRNA.paralogs_removed.bam
	""" > ../bsubFiles/flu_capseq_bowtie_${time}.paralogs_removed.bsub
done
