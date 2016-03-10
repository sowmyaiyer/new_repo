# ./processChipNew_pe_JK_vs_DC_preprocess.sh
#for sample in  `awk '{ print $NF}' ../txt/chip01132016_sampleInfo.txt`
for sample in  `awk '{ print $NF}' ../txt/chip01132016_sampleInfo_H3K27ac_only.txt`
do
	echo $sample
	expression_for_celltype=`echo $sample | awk '{split($1,arr,"_"); if (substr(arr[4],1,2) == "JK") print "gene_tss_expression_Jurkat"; else print "gene_tss_expression_time_0h"}'`
	echo """
	module load IGVTools/2.3.31
	r1_lines=\`zcat /farline/umw_garberlab/human/DC/ChIP-Seq/ChIP_seq_1_13_2016/${sample}.R1.fastq.gz | wc -l\`
	r2_lines=\`zcat /farline/umw_garberlab/human/DC/ChIP-Seq/ChIP_seq_1_13_2016/${sample}.R2.fastq.gz | wc -l\`
	if [[ \$r1_lines != \$r2_lines ]]; then
		echo ${sample} has unequal reads for R1 and R2
		exit 1
	fi
	bwa aln -t 8 /share/data/umw_biocore/Genomes/human/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.6.0/genome.fa  /farline/umw_garberlab/human/DC/ChIP-Seq/ChIP_seq_1_13_2016/${sample}.R1.fastq.gz > ../bwa_out/${sample}.pe.R1.sai
	bwa aln -t 8 /share/data/umw_biocore/Genomes/human/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.6.0/genome.fa  /farline/umw_garberlab/human/DC/ChIP-Seq/ChIP_seq_1_13_2016/${sample}.R2.fastq.gz > ../bwa_out/${sample}.pe.R2.sai
	bwa sampe /share/data/umw_biocore/Genomes/human/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.6.0/genome.fa  ../bwa_out/${sample}.pe.R1.sai  ../bwa_out/${sample}.pe.R2.sai /farline/umw_garberlab/human/DC/ChIP-Seq/ChIP_seq_1_13_2016/${sample}.R1.fastq.gz /farline/umw_garberlab/human/DC/ChIP-Seq/ChIP_seq_1_13_2016/${sample}.R2.fastq.gz > ../bwa_out/${sample}.pe.sam
	samtools view -bS -F 4 -o ../bwa_out/${sample}.pe.bam ../bwa_out/${sample}.pe.sam
	samtools view -bS -f 2 -o ../bwa_out/${sample}.pe.properly_paired.bam ../bwa_out/${sample}.pe.sam
	samtools sort -T ../bwa_out/${sample}.pe.sorted -o ../bwa_out/${sample}.pe.sorted.bam ../bwa_out/${sample}.pe.bam
	samtools sort -T ../bwa_out/${sample}.pe.properly_paired.sorted -o ../bwa_out/${sample}.pe.properly_paired.sorted.bam ../bwa_out/${sample}.pe.properly_paired.bam
        samtools index ../bwa_out/${sample}.pe.sorted.bam
        samtools index ../bwa_out/${sample}.pe.properly_paired.sorted.bam
	rm   ../bwa_out/${sample}.pe.R1.sai ../bwa_out/${sample}.pe.R2.sai ../bwa_out/${sample}.pe.sam ../bwa_out/${sample}.pe.bam ../bwa_out/${sample}.pe.properly_paired.bam
	java -jar /share/pkg/picard/1.96/DownsampleSam.jar INPUT=../bwa_out/${sample}.pe.properly_paired.sorted.bam PROBABILITY=0.25 RANDOM_SEED=null OUTPUT=../bwa_out/${sample}.pe.properly_paired.downsampled.bam VALIDATION_STRINGENCY=LENIENT
        samtools sort ../bwa_out/${sample}.pe.properly_paired.downsampled.bam -o ../bwa_out/${sample}.pe.properly_paired.downsampled.sorted.bam -T ../bwa_out/${sample}.pe.properly_paired.downsampled.sorted
#        # begin QC
        samtools flagstat ../bwa_out/${sample}.pe.properly_paired.downsampled.sorted.bam > ../bwa_out/${sample}.downsampled.flagstat
        mappedReads=\`samtools view  -c ../bwa_out/${sample}.pe.properly_paired.downsampled.sorted.bam\`
        mappedReadsInPromoter=\`samtools view -F 4 -b ../bwa_out/${sample}.pe.properly_paired.downsampled.sorted.bam | bedtools bamtobed -i stdin | bedtools intersect -a stdin -b ../txt/hg19_promoter.bed -u | wc -l\`
        echo mapped reads in promoter = \$mappedReadsInPromoter >> ../bwa_out/${sample}.downsampled.flagstat
	num_properly_paired=\`grep \"properly paired\" ../bwa_out/${sample}.downsampled.flagstat | cut -d\" \" -f1\`
	java -jar /share/pkg/picard/1.96/CollectInsertSizeMetrics.jar INPUT=../bwa_out/${sample}.pe.properly_paired.downsampled.sorted.bam HISTOGRAM_FILE=../results/${sample}.downsampled.pdf OUTPUT=../results/${sample}.downsampled.txt
	java -jar /share/pkg/picard/1.96/MarkDuplicates.jar INPUT=../bwa_out/${sample}.pe.properly_paired.downsampled.sorted.bam METRICS_FILE=../results/${sample}.pcrdups.properly_paired.downsampled.txt REMOVE_DUPLICATES=true OUTPUT=../results/delete.${sample}.dups_removed.properly_paired.downsampled.bam
	rm ../results/delete.${sample}.dups_removed.properly_paired.downsampled.bam 
	percent_dups=\`grep \"Unknown Library\" ../results/${sample}.pcrdups.properly_paired.downsampled.txt | awk '{ print \$(NF-1)}'\`
	totalreads=\$((\${r1_lines} / 4))
	echo ${sample},\${totalreads}, \${mappedReads}, \${num_properly_paired},\${mappedReadsInPromoter},\${percent_dups} > ../txt/${sample}.downsampled.stats
        scalingFactor=\`echo \$mappedReads | awk -vmappedReads=\$mappedReads 'BEGIN{print 1000000/mappedReads}'\`
        bedtools genomecov -ibam ../bwa_out/${sample}.pe.properly_paired.downsampled.sorted.bam -g ~/gnearline/common_scripts/hg19.genome -bga -scale \${scalingFactor} > ../bwa_out/${sample}.downsampled.bedGraph
        bedGraphToBigWig ../bwa_out/${sample}.downsampled.bedGraph ~/gnearline/common_scripts/hg19.sorted.genome ../bwa_out/${sample}.downsampled.bw
	rm ../bwa_out/${sample}.downsampled.bedGraph
        for exp in {\"high\",\"med\",\"low\",\"verylow\"}
        do
                bedtools slop -i ../txt/${expression_for_celltype}.\${exp}.minus.bed -b 1000 -g ~/gnearline/common_scripts/hg19.genome | bedtools makewindows -b stdin -w 10 -reverse -i srcwinnum | bigWigAverageOverBed ../bwa_out/${sample}.downsampled.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/${sample}.\${exp}.downsampled.minus.bed
                bedtools slop -i ../txt/${expression_for_celltype}.\${exp}.plus.bed -b 1000 -g ~/gnearline/common_scripts/hg19.genome | bedtools makewindows -b stdin -w 10 -i srcwinnum | bigWigAverageOverBed ../bwa_out/${sample}.downsampled.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/${sample}.\${exp}.downsampled.plus.bed
                sed 's/_/\t/g' ../txt/${sample}.\${exp}.downsampled.minus.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' > ../txt/${sample}.\${exp}.downsampled.all_signals.txt
                sed 's/_/\t/g' ../txt/${sample}.\${exp}.downsampled.plus.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' >> ../txt/${sample}.\${exp}.downsampled.all_signals.txt
        done
        Rscript aggrPlots.R ../txt/${sample}.high.downsampled.all_signals.txt ../txt/${sample}.med.downsampled.all_signals.txt ../txt/${sample}.low.downsampled.all_signals.txt ../txt/${sample}.verylow.downsampled.all_signals.txt ../results/${sample}.downsampled.aggr.pdf $sample
	igvtools count -z 5 -w 25 -e 250 ../bwa_out/${sample}.pe.properly_paired.downsampled.sorted.bam  ../bwa_out/${sample}.pe.properly_paired.downsampled.sorted.tdf  hg19
	""" > ../bsubFiles/new_JK_vs_DC_pe_chip01132016_${sample}.downsampled.bsub
done
