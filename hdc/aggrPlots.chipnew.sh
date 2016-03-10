time="0h"
for sample in  `awk '{ print $1}' ../txt/chipnew_sampleInfo.txt | sort | uniq`
do
	echo """
	samtools flagstat ../bwa_out/${sample}.se.sorted.bam > ../bwa_out/${sample}.flagstat
	mappedReads=\`samtools view -F 0x4  ../bwa_out/${sample}.se.sorted.bam | wc -l\`
	mappedReadsInPromoter=\`samtools view -F 0x4 ../bwa_out/${sample}.se.sorted.bam | bedtools bamtobed -i stdin | bedtools intersect -a stdin -b ../txt/hg19_promoter.bed -u | wc -l\`
	scalingFactor=\`echo \$mappedReads | awk -vmappedReads=\$mappedReads 'BEGIN{print 1000000/mappedReads}'\`
	bedtools genomecov -ibam ../bwa_out/${sample}.se.sorted.bam -g ~/gnearline/common_scripts/hg19.genome -bga -scale \${scalingFactor} > ../bwa_out/${sample}.bedGraph
	bedGraphToBigWig ../bwa_out/${sample}.bedGraph ~/gnearline/common_scripts/hg19.genome ../bwa_out/${sample}.bw
	for exp in {\"high\",\"med\",\"low\",\"verylow\"}
	do
		bedtools slop -i ../txt/gene_tss_expression_time_${time}.\${exp}.minus.bed -b 1000 -g ~/gnearline/common_scripts/hg19.genome | bedtools makewindows -b stdin -w 10 -reverse -i srcwinnum | bigWigAverageOverBed ../bwa_out/${sample}.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/${sample}.\${exp}.minus.bed
		bedtools slop -i ../txt/gene_tss_expression_time_${time}.\${exp}.plus.bed -b 1000 -g ~/gnearline/common_scripts/hg19.genome | bedtools makewindows -b stdin -w 10 -i srcwinnum | bigWigAverageOverBed ../bwa_out/${sample}.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/${sample}.\${exp}.plus.bed
		sed 's/_/\t/g' ../txt/${sample}.\${exp}.minus.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' > ../txt/${sample}.\${exp}.all_signals.txt	
		sed 's/_/\t/g' ../txt/${sample}.\${exp}.plus.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' >> ../txt/${sample}.\${exp}.all_signals.txt
	done
	Rscript aggrPlots.R ../txt/${sample}.high.all_signals.txt ../txt/${sample}.med.all_signals.txt ../txt/${sample}.low.all_signals.txt ../txt/${sample}.verylow.all_signals.txt ../results/${sample}.H3K4me3.pdf 
	""" > ../bsubFiles/qc_chipnew_${sample}.bsub
done
