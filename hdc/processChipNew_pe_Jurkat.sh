#while read line
#do
#        r1_file=`echo $line | awk '{ print $2}'`
#        r2_file=`echo $line | awk '{ print $3}'`
#        sample=`echo $line | awk '{ print $1}'`
 #       cat /farline/umw_garberlab/mouse/DC/ATAC-Seq/ATAC_CHIP_12_15_15/fastq.untrimmed/${r1_file} >> ~/gfarline/hdc/Jurkat_new/${sample}_R1.gz
#        cat /farline/umw_garberlab/mouse/DC/ATAC-Seq/ATAC_CHIP_12_15_15/fastq.untrimmed/${r2_file} >> ~/gfarline/hdc/Jurkat_new/${sample}_R2.gz
#        ls -ltr ~/gfarline/hdc/Jurkat_new/${sample}_R[1,2].gz
#done<../txt/chip1215_sampleInfo.txt

for sample in  `awk '{ print $1}' ../txt/chip1215_sampleInfo.txt | sort | uniq`
do
	echo $sample
	echo """
	r1_lines=\`zcat ~/gfarline/hdc/Jurkat_new/${sample}_R1.gz | wc -l\`
	r2_lines=\`zcat ~/gfarline/hdc/Jurkat_new/${sample}_R2.gz | wc -l\`
	if [[ \$r1_lines != \$r2_lines ]]; then
		echo ${sample} has unequal reads for R1 and R2
		exit 1
	fi
	bwa aln -t 8 /share/data/umw_biocore/Genomes/human/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.6.0/genome.fa  ~/gfarline/hdc/Jurkat_new/${sample}_R1.gz > ../bwa_out/${sample}.pe.R1.sai
	bwa aln -t 8 /share/data/umw_biocore/Genomes/human/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.6.0/genome.fa  ~/gfarline/hdc/Jurkat_new/${sample}_R2.gz > ../bwa_out/${sample}.pe.R2.sai
	bwa sampe /share/data/umw_biocore/Genomes/human/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.6.0/genome.fa  ../bwa_out/${sample}.pe.R1.sai  ../bwa_out/${sample}.pe.R2.sai ~/gfarline/hdc/Jurkat_new/${sample}_R1.gz ~/gfarline/hdc/Jurkat_new/${sample}_R2.gz > ../bwa_out/${sample}.pe.sam
	samtools view -bS -F 4 -o ../bwa_out/${sample}.pe.bam ../bwa_out/${sample}.pe.sam
	samtools view -bS -f 2 -o ../bwa_out/${sample}.pe.properly_paired.bam ../bwa_out/${sample}.pe.sam
	samtools sort -T ../bwa_out/${sample}.pe.sorted -o ../bwa_out/${sample}.pe.sorted.bam ../bwa_out/${sample}.pe.bam
	samtools sort -T ../bwa_out/${sample}.pe.properly_paired.sorted -o ../bwa_out/${sample}.pe.properly_paired.sorted.bam ../bwa_out/${sample}.pe.properly_paired.bam
        samtools index ../bwa_out/${sample}.pe.sorted.bam
        samtools index ../bwa_out/${sample}.pe.properly_paired.sorted.bam
	rm   ../bwa_out/${sample}.pe.R1.sai ../bwa_out/${sample}.pe.R2.sai ../bwa_out/${sample}.pe.sam ../bwa_out/${sample}.pe.bam ../bwa_out/${sample}.pe.properly_paired.bam

#        # begin QC
        samtools flagstat ../bwa_out/${sample}.pe.sorted.bam > ../bwa_out/${sample}.flagstat
        mappedReads=\`samtools view  ../bwa_out/${sample}.pe.properly_paired.sorted.bam | wc -l\`
        mappedReadsInPromoter=\`samtools view -F 4 -b ../bwa_out/${sample}.pe.properly_paired.sorted.bam | bedtools bamtobed -i stdin | bedtools intersect -a stdin -b ../txt/hg19_promoter.bed -u | wc -l\`
        echo mapped reads in promoter = \$mappedReadsInPromoter >> ../bwa_out/${sample}.flagstat
	
	java -jar /share/pkg/picard/1.96/CollectInsertSizeMetrics.jar INPUT=/home/si14w/gnearline/hdc/bwa_out/${sample}.pe.sorted.bam HISTOGRAM_FILE=../results/${sample}.pdf OUTPUT=../results/${sample}.txt
	java -jar /share/pkg/picard/1.96/MarkDuplicates.jar INPUT=/home/si14w/gnearline/hdc/bwa_out/${sample}.pe.properly_paired.sorted.bam METRICS_FILE=../results/${sample}.pcrdups.properly_paired.txt REMOVE_DUPLICATES=true OUTPUT=../results/delete.${sample}.dups_removed.properly_paired.bam
	java -jar /share/pkg/picard/1.96/MarkDuplicates.jar INPUT=/home/si14w/gnearline/hdc/bwa_out/${sample}.pe.sorted.bam METRICS_FILE=../results/${sample}.pcrdups.txt REMOVE_DUPLICATES=true OUTPUT=../results/delete.${sample}.dups_removed.bam
	rm ../results/delete.${sample}.dups_removed.bam ../results/delete.${sample}.dups_removed.properly_paired.bam
	
        scalingFactor=\`echo \$mappedReads | awk -vmappedReads=\$mappedReads 'BEGIN{print 1000000/mappedReads}'\`
        bedtools genomecov -ibam ../bwa_out/${sample}.pe.properly_paired.sorted.bam -g ~/gnearline/common_scripts/hg19.genome -bga -scale \${scalingFactor} > ../bwa_out/${sample}.bedGraph
        bedGraphToBigWig ../bwa_out/${sample}.bedGraph ~/gnearline/common_scripts/hg19.sorted.genome ../bwa_out/${sample}.bw
	rm ../bwa_out/${sample}.bedGraph
        for exp in {\"high\",\"med\",\"low\",\"verylow\"}
        do
                bedtools slop -i ../txt/gene_tss_expression_Jurkat.\${exp}.minus.bed -b 1000 -g ~/gnearline/common_scripts/hg19.genome | bedtools makewindows -b stdin -w 10 -reverse -i srcwinnum | bigWigAverageOverBed ../bwa_out/${sample}.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/${sample}.\${exp}.minus.bed
                bedtools slop -i ../txt/gene_tss_expression_Jurkat.\${exp}.plus.bed -b 1000 -g ~/gnearline/common_scripts/hg19.genome | bedtools makewindows -b stdin -w 10 -i srcwinnum | bigWigAverageOverBed ../bwa_out/${sample}.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/${sample}.\${exp}.plus.bed
                sed 's/_/\t/g' ../txt/${sample}.\${exp}.minus.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' > ../txt/${sample}.\${exp}.all_signals.txt
                sed 's/_/\t/g' ../txt/${sample}.\${exp}.plus.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' >> ../txt/${sample}.\${exp}.all_signals.txt
        done
        Rscript aggrPlots.R ../txt/${sample}.high.all_signals.txt ../txt/${sample}.med.all_signals.txt ../txt/${sample}.low.all_signals.txt ../txt/${sample}.verylow.all_signals.txt ../results/${sample}.aggr.pdf
	igvtools count -z 5 -w 25 -e 250 ../bwa_out/${sample}.pe.properly_paired.sorted.bam  ../bwa_out/${sample}.pe.properly_paired.tdf  hg19
	""" > ../bsubFiles/Jurkat_pe_chip1215_${sample}.bsub
done
