time="0h"
for sample in  `awk '{ print $1}' ../txt/chip1124_sampleInfo.txt | sort | uniq`
do
	echo $sample
	echo """
	r1_lines=\`zcat ~/gnearline/hdc/chip1124/${sample}_R1.gz | wc -l\`
	r2_lines=\`zcat ~/gnearline/hdc/chip1124/${sample}_R2.gz | wc -l\`
	if [[ \$r1_lines != \$r2_lines ]]; then
		echo ${sample} has unequal reads for R1 and R2
		exit 1
	fi
	bwa aln -t 8 /share/data/umw_biocore/Genomes/human/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.6.0/genome.fa  ~/gnearline/hdc/chip1124/${sample}_R1.gz > ../bwa_out/${sample}.pe.R1.sai
	bwa aln -t 8 /share/data/umw_biocore/Genomes/human/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.6.0/genome.fa  ~/gnearline/hdc/chip1124/${sample}_R2.gz > ../bwa_out/${sample}.pe.R2.sai
	bwa sampe /share/data/umw_biocore/Genomes/human/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.6.0/genome.fa  ../bwa_out/${sample}.pe.R1.sai  ../bwa_out/${sample}.pe.R2.sai ~/gnearline/hdc/chip1124/${sample}_R1.gz ~/gnearline/hdc/chip1124/${sample}_R2.gz > ../bwa_out/${sample}.pe.sam
	samtools view -bS -F 4 -o ../bwa_out/${sample}.pe.bam ../bwa_out/${sample}.pe.sam
	samtools sort -T ../bwa_out/${sample}.pe.sorted -o ../bwa_out/${sample}.pe.sorted.bam ../bwa_out/${sample}.pe.bam
        samtools index ../bwa_out/${sample}.pe.sorted.bam
	rm   ../bwa_out/${sample}.pe.R1.sai ../bwa_out/${sample}.pe.R2.sai ../bwa_out/${sample}.pe.sam

#        # begin QC
        samtools flagstat ../bwa_out/${sample}.pe.sorted.bam > ../bwa_out/${sample}.flagstat
        mappedReads=\`samtools view -F 4  ../bwa_out/${sample}.pe.sorted.bam | wc -l\`
        mappedReadsInPromoter=\`samtools view -F 4 -b ../bwa_out/${sample}.pe.sorted.bam | bedtools bamtobed -i stdin | bedtools intersect -a stdin -b ../txt/hg19_promoter.bed -u | wc -l\`
        echo mapped reads in promoter = \$mappedReadsInPromoter >> ../bwa_out/${sample}.flagstat
	
	java -jar /share/pkg/picard/1.96/CollectInsertSizeMetrics.jar INPUT=/home/si14w/gnearline/hdc/bwa_out/${sample}.pe.sorted.bam HISTOGRAM_FILE=../results/${sample}.${time}.pdf OUTPUT=../results/${sample}.${time}.txt
	java -jar /share/pkg/picard/1.96/MarkDuplicates.jar INPUT=/home/si14w/gnearline/hdc/bwa_out/${sample}.pe.sorted.bam METRICS_FILE=../results/${sample}.${time}.pcrdups.txt REMOVE_DUPLICATES=true
	
        scalingFactor=\`echo \$mappedReads | awk -vmappedReads=\$mappedReads 'BEGIN{print 1000000/mappedReads}'\`
        bedtools genomecov -ibam ../bwa_out/${sample}.pe.sorted.bam -g ~/gnearline/common_scripts/hg19.genome -bga -scale \${scalingFactor} > ../bwa_out/${sample}.bedGraph
        bedGraphToBigWig ../bwa_out/${sample}.bedGraph ~/gnearline/common_scripts/hg19.genome ../bwa_out/${sample}.bw
        for exp in {\"high\",\"med\",\"low\",\"verylow\"}
        do
                bedtools slop -i ../txt/gene_tss_expression_time_${time}.\${exp}.minus.bed -b 1000 -g ~/gnearline/common_scripts/hg19.genome | bedtools makewindows -b stdin -w 10 -reverse -i srcwinnum | bigWigAverageOverBed ../bwa_out/${sample}.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/${sample}.\${exp}.minus.bed
                bedtools slop -i ../txt/gene_tss_expression_time_${time}.\${exp}.plus.bed -b 1000 -g ~/gnearline/common_scripts/hg19.genome | bedtools makewindows -b stdin -w 10 -i srcwinnum | bigWigAverageOverBed ../bwa_out/${sample}.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/${sample}.\${exp}.plus.bed
                sed 's/_/\t/g' ../txt/${sample}.\${exp}.minus.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' > ../txt/${sample}.\${exp}.all_signals.txt
                sed 's/_/\t/g' ../txt/${sample}.\${exp}.plus.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' >> ../txt/${sample}.\${exp}.all_signals.txt
        done
	sed -i '/GTF2H2C/d' ../txt/${sample}.med.all_signals.txt
        Rscript aggrPlots.R ../txt/${sample}.high.all_signals.txt ../txt/${sample}.med.all_signals.txt ../txt/${sample}.low.all_signals.txt ../txt/${sample}.verylow.all_signals.txt ../results/${sample}.H3K4me3.pdf
	""" > ../bsubFiles/bwa_pe_chip1124_${sample}.bsub
done
