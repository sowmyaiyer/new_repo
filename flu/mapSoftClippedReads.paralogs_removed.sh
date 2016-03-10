module load bowtie2/2-2.1.0
awk '{ if ($6 == "-") {tss=$3} else if ($6 == "+") {tss=$2} printf("%s\t%d\t%d\t%s\t%s\t%s\n", $1,tss,tss,$4,$5,$6) }' ../txt/hg19_one_tr_per_gene.bed | bedtools slop -i stdin -g ~/gnearline/common_scripts/hg19.genome -b 50 | bedtools getfasta -name -fi /share/data/umw_biocore/genome_data/human/hg19full/hg19full.fa -bed stdin -fo ../txt/hg19_one_tr_per_gene.fa
awk '{ if ($6 == "-") { tss=$3} else if ($6 == "+") { tss=$2} printf("%s\t%d\t%d\t%s_%s_%d_%d\t%s\t%s\t%s\n",$1,tss,tss,$4,$1,$2,$3,$5,$6,$7)}' ../txt/hg19_snRNA_and_snoRNA.tmp.bed | bedtools slop -i stdin -g ~/gnearline/common_scripts/hg19.genome -b 50 | bedtools getfasta -name -fi /share/data/umw_biocore/genome_data/human/hg19full/hg19full.fa -bed stdin -fo ../txt/hg19_snRNA_and_snoRNA.sequences.fa
./removesnRNAParalogs.sh
cat ../txt/hg19_snRNA_and_snoRNA.sequences.paralogs_removed.fa ../txt/hg19_one_tr_per_gene.fa > ../txt/hg19_one_tr_per_gene_and_snRNA_snoRNA.paralogs_removed.fa
bowtie2-build ../txt/hg19_one_tr_per_gene_and_snRNA_snoRNA.paralogs_removed.fa ../bowtie_out/hg19_one_tr_per_gene_and_snRNA_snoRNA.paralogs_removed.fa
for time in {"A","B","C","D","E","F","G","H"}
do
	samtools view ../bowtie_out/112213flucap_${time}.sorted.bam | grep "B59FULL" | awk '
	{ 
		if ($4 == 1) 
		{
			split($6,arr,"S"); 
			if(arr[1] >= 9 && arr[1] <= 15) 
			printf(">%s\n%s\n",$1,substr($10,1,arr[1]))
		} 
	}' > ../bowtie_out/112213flucap_${time}.flu_softclipped.fa
	echo """
	module load bowtie2/2-2.1.0
	bowtie2 -p 8 --met-file ../bowtie_out/softClip_map_stats.${time}.paralogs_removed.txt -x ../bowtie_out/hg19_one_tr_per_gene_and_snRNA_snoRNA.paralogs_removed.fa -S --norc -N 0 -U ../bowtie_out/112213flucap_${time}.flu_softclipped.fa --all -f -S ../bowtie_out/112213flucap_softclip.${time}.paralogs_removed.sam
        samtools view -bS -F 4 ../bowtie_out/112213flucap_softclip.${time}.paralogs_removed.sam > ../bowtie_out/112213flucap_softclip.${time}.paralogs_removed.bam
        samtools view -b -F 16 ../bowtie_out/112213flucap_softclip.${time}.paralogs_removed.bam | samtools sort -n -T ../bowtie_out/112213flucap_${time}.paralogs_removed.sorted -o ../bowtie_out/112213flucap_${time}.softclip.sorted.paralogs_removed.bam -
	samtools view ../bowtie_out/112213flucap_${time}.softclip.sorted.paralogs_removed.bam | grep "XM:i:0" | awk '{ print \$1}' | sort | uniq -c | awk '{ print \$1}' >  ../txt/histogram_multireads_${time}.paralogs_removed.txt
        rm ../bowtie_out/112213flucap_softclip.${time}.paralogs_removed.sam ../bowtie_out/112213flucap_softclip.${time}.paralogs_removed.bam
	""" > ../bsubFiles/mapAndProcess_capseq_softclip.${time}.paralogs_removed.bsub
#        bedtools intersect -nonamecheck -a ../txt/hg19_genes_150.bed -b ../bowtie_out/112213flucap_softclip.${time}.bed -c -F 1 > ../txt/112213flucap_${time}.readCounts_150bp.allRNA.softclipped.genes.bed
#        Rscript normalizeReadCounts.R $HOME/gnearline/flu/txt/112213flucap_${time}.readCounts_150bp.allRNA.softclipped.genes.bed $HOME/gnearline/flu/txt/112213flucap_${time}.readCounts_150bp.rpm.allRNA.softclipped.genes.bed
done
