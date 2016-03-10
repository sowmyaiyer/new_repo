#module load fastqc/0.10.1
#while read line
#do
#        cat /farline/umw_garberlab/human/DC/ChIP-Seq/ChIP_seq_12_22_2015/${line}*R1_001.fastq.gz > ~/gfarline/hdc/DC_new/${line}_R1.gz
#        cat /farline/umw_garberlab/human/DC/ChIP-Seq/ChIP_seq_12_22_2015/${line}*R2_001.fastq.gz > ~/gfarline/hdc/DC_new/${line}_R2.gz
#        ls -ltr ~/gfarline/hdc/DC_new/${line}_R[1,2].gz
#done<../txt/chip1222_sampleInfo_uniq.txt
#for  fastq_gz in  /home/si14w/gfarline/hdc/DC_new/*gz
#do
#         fastqc $fastq_gz -o /home/si14w/gnearline/hdc/chipseq_fastqc
#done

for sample in  `awk '{ print $1}' ../txt/chip1222_sampleInfo_uniq.txt`
do
	echo $sample
	echo """
	r1_lines=\`zcat ~/gfarline/hdc/DC_new/${sample}_R1.gz | wc -l\`
	r2_lines=\`zcat ~/gfarline/hdc/DC_new/${sample}_R2.gz | wc -l\`
	if [[ \$r1_lines != \$r2_lines ]]; then
		echo ${sample} has unequal reads for R1 and R2
		exit 1
	fi
	bwa aln -t 8 /share/data/umw_biocore/Genomes/human/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.6.0/genome.fa  ~/gfarline/hdc/DC_new/${sample}_R1.gz > ../bwa_out/${sample}.pe.R1.sai
	bwa aln -t 8 /share/data/umw_biocore/Genomes/human/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.6.0/genome.fa  ~/gfarline/hdc/DC_new/${sample}_R2.gz > ../bwa_out/${sample}.pe.R2.sai
	bwa sampe /share/data/umw_biocore/Genomes/human/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.6.0/genome.fa  ../bwa_out/${sample}.pe.R1.sai  ../bwa_out/${sample}.pe.R2.sai ~/gfarline/hdc/DC_new/${sample}_R1.gz ~/gfarline/hdc/DC_new/${sample}_R2.gz > ../bwa_out/${sample}.pe.sam
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
                bedtools slop -i ../txt/gene_tss_expression_time_0h.\${exp}.minus.bed -b 1000 -g ~/gnearline/common_scripts/hg19.genome | bedtools makewindows -b stdin -w 10 -reverse -i srcwinnum | bigWigAverageOverBed ../bwa_out/${sample}.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/${sample}.\${exp}.minus.bed
                bedtools slop -i ../txt/gene_tss_expression_time_0h.\${exp}.plus.bed -b 1000 -g ~/gnearline/common_scripts/hg19.genome | bedtools makewindows -b stdin -w 10 -i srcwinnum | bigWigAverageOverBed ../bwa_out/${sample}.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/${sample}.\${exp}.plus.bed
                sed 's/_/\t/g' ../txt/${sample}.\${exp}.minus.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' > ../txt/${sample}.\${exp}.all_signals.txt
                sed 's/_/\t/g' ../txt/${sample}.\${exp}.plus.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' >> ../txt/${sample}.\${exp}.all_signals.txt
        done
        Rscript aggrPlots.R ../txt/${sample}.high.all_signals.txt ../txt/${sample}.med.all_signals.txt ../txt/${sample}.low.all_signals.txt ../txt/${sample}.verylow.all_signals.txt ../results/${sample}.aggr.pdf
	igvtools count -z 5 -w 25 -e 250 ../bwa_out/${sample}.pe.properly_paired.sorted.bam  ../bwa_out/${sample}.pe.properly_paired.tdf  hg19
	""" > ../bsubFiles/new_DC_pe_chip1222_${sample}.bsub
done
