for ab in {"H3K4me3","H3K27ac","H3K4me1","IgG","input"}
do
	for time in {"0h","2h_LPS","4h_LPS"}
	do
		echo """
#		samtools flagstat /project/umw_garberlab/human/DC/ChIP-Seq/ChIP_seq_08_05_2015/pipeline/seqmapping/chip/E01_6d_${time}_${ab}.sorted.bam > ../txt/E01_6d_${time}_${ab}.flagstat
		reads_mapped=\`awk '{ if (NR == 5) print \$1}' ../txt/E01_6d_${time}_${ab}.flagstat\`
		reads_properly_paired=\`awk '{ if (NR == 9) print \$1}' ../txt/E01_6d_${time}_${ab}.flagstat\`
		reads_mate_diff_chr=\`awk '{ if (NR == 12) print \$1}' ../txt/E01_6d_${time}_${ab}.flagstat\`
#		bedtools bamtobed -i /project/umw_garberlab/human/DC/ChIP-Seq/ChIP_seq_08_05_2015/pipeline/seqmapping/chip/E01_6d_${time}_${ab}.sorted.bam > ../txt/E01_6d_${time}_${ab}.sorted.bam.bamtobed.bed
		reads_in_promoter=\`bedtools slop -i /home/si14w/nearline/enhancer_predictions/gencodeV19_TSS.bed -b 1000 -g ~/gnearline/common_scripts/hg19.genome | bedtools intersect -a ../txt/E01_6d_${time}_${ab}.sorted.bam.bamtobed.bed -b stdin -u | wc -l\`
		num_peaks=\`wc -l /project/umw_garberlab/human/DC/ChIP-Seq/ChIP_seq_08_05_2015/pipeline/macs/E01_6d_${time}_${ab}.vs.E01_6d_${time}_input_peaks.bed \`
		frip=\`bedtools intersect -a ../txt/E01_6d_${time}_${ab}.sorted.bam.bamtobed.bed -b /project/umw_garberlab/human/DC/ChIP-Seq/ChIP_seq_08_05_2015/pipeline/macs/E01_6d_${time}_${ab}.vs.E01_6d_${time}_input_peaks.bed -u | wc -l\`
		samtools view -f 0x02  /project/umw_garberlab/human/DC/ChIP-Seq/ChIP_seq_08_05_2015/pipeline/seqmapping/chip/E01_6d_${time}_${ab}.sorted.bam -b | samtools sort -o ../txt/E01_6d_${time}_${ab}.properly_paired.sorted.bam -T ../txt/E01_6d_${time}_${ab}.sorted
		java -jar /share/pkg/picard/1.96/CollectInsertSizeMetrics.jar INPUT=../txt/E01_6d_${time}_${ab}.properly_paired.sorted.bam HISTOGRAM_FILE=../results/qc_insert_size_${time}_${ab}.properly_paired.pdf OUTPUT=../results/qc_insert_size_${time}_${ab}.properly_paired.txt
		java -jar /share/pkg/picard/1.96/CollectInsertSizeMetrics.jar INPUT=/project/umw_garberlab/human/DC/ChIP-Seq/ChIP_seq_08_05_2015/pipeline/seqmapping/chip/E01_6d_${time}_${ab}.sorted.bam  HISTOGRAM_FILE=../results/qc_insert_size_${time}_${ab}.allreads.pdf OUTPUT=../results/qc_insert_size_${time}_${ab}.allreads.txt
		bedtools bamtobed -i ../txt/E01_6d_${time}_${ab}.properly_paired.sorted.bam > ../txt/E01_6d_${time}_${ab}.properly_paired.sorted.bamtobed.bed
		properly_paired_reads_in_promoter=\`bedtools slop -i /home/si14w/nearline/enhancer_predictions/gencodeV19_TSS.bed -b 1000 -g ~/gnearline/common_scripts/hg19.genome | bedtools intersect -a ../txt/E01_6d_${time}_${ab}.properly_paired.sorted.bamtobed.bed  -b stdin -u | wc -l\`
		echo ${ab} ${time} \${reads_mapped} \${reads_properly_paired} \${reads_mate_diff_chr} \${reads_in_promoter} \${properly_paired_reads_in_promoter} \${num_peaks} \${frip} > ../txt/qc_${time}_${ab}.txt
		""" > ../bsubFiles/qc_${time}_${ab}.bsub
	done
	echo antibody time reads_mapped reads_properly_paired reads_mate_diff_chr reads_in_promoter properly_paired_reads_in_promoter num_peaks frip > ../txt/qc_all.txt
	cat ../txt/qc_*h*txt >> ../txt/qc_all.txt
done

