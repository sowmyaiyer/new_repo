# ./getGeneExpressionData.sh
# ./processChipNew_pe_JK_vs_DC_preprocess.sh

while read sample
do
#	time_tmp=`echo $sample | cut -d"_" -f2`
#        time_tmp2=${time_tmp}"_LPS"
#        if [[ ($time_tmp == "CV") || ($time_tmp == "PM") || ($time_tmp == "30m") ]]; then
                time="0h"
#        else
#                time=$time_tmp2
#        fi
#        echo $time_tmp2 $time
	echo """
	module load cutadapt/1.9
	module load IGVTools/2.3.31

	cutadapt --minimum-length 30 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o /farline/umw_garberlab/human/DC/ChIP-Seq/patient_samples/chip_retest_2/${sample}.R1.trimmed.fastq.gz -p  /farline/umw_garberlab/human/DC/ChIP-Seq/patient_samples/chip_retest_2/${sample}.R2.trimmed.fastq.gz  /farline/umw_garberlab/human/DC/ChIP-Seq/patient_samples/chip_retest_2/${sample}.R1.fastq.gz  /farline/umw_garberlab/human/DC/ChIP-Seq/patient_samples/chip_retest_2/${sample}.R2.fastq.gz

	r1_lines=\`zcat /farline/umw_garberlab/human/DC/ChIP-Seq/patient_samples/chip_retest_2/${sample}.R1.trimmed.fastq.gz | wc -l\`
	r2_lines=\`zcat /farline/umw_garberlab/human/DC/ChIP-Seq/patient_samples/chip_retest_2/${sample}.R2.trimmed.fastq.gz | wc -l\`
	if [[ \$r1_lines != \$r2_lines ]]; then
		echo ${sample} has unequal reads for R1 and R2
		exit 1
	fi
	bwa aln -t 8 /share/data/umw_biocore/Genomes/human/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.6.0/genome.fa  /farline/umw_garberlab/human/DC/ChIP-Seq/patient_samples/chip_retest_2/${sample}.R1.trimmed.fastq.gz > ../bwa_out/${sample}.pe.R1.sai
	bwa aln -t 8 /share/data/umw_biocore/Genomes/human/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.6.0/genome.fa  /farline/umw_garberlab/human/DC/ChIP-Seq/patient_samples/chip_retest_2/${sample}.R2.trimmed.fastq.gz > ../bwa_out/${sample}.pe.R2.sai
	bwa sampe /share/data/umw_biocore/Genomes/human/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.6.0/genome.fa  ../bwa_out/${sample}.pe.R1.sai  ../bwa_out/${sample}.pe.R2.sai /farline/umw_garberlab/human/DC/ChIP-Seq/patient_samples/chip_retest_2/${sample}.R1.trimmed.fastq.gz /farline/umw_garberlab/human/DC/ChIP-Seq/patient_samples/chip_retest_2/${sample}.R2.trimmed.fastq.gz > ../bwa_out/${sample}.pe.sam
	samtools view -bS -F 4 -o ../bwa_out/${sample}.pe.bam ../bwa_out/${sample}.pe.sam
	samtools view -bS -f 2 -o ../bwa_out/${sample}.pe.properly_paired.bam ../bwa_out/${sample}.pe.sam
	samtools sort -T ../bwa_out/${sample}.pe.sorted -o ../bwa_out/${sample}.pe.sorted.bam ../bwa_out/${sample}.pe.bam
	samtools sort -T ../bwa_out/${sample}.pe.properly_paired.sorted -o ../bwa_out/${sample}.pe.properly_paired.sorted.bam ../bwa_out/${sample}.pe.properly_paired.bam
        samtools index ../bwa_out/${sample}.pe.sorted.bam
        samtools index ../bwa_out/${sample}.pe.properly_paired.sorted.bam
	rm   ../bwa_out/${sample}.pe.R1.sai ../bwa_out/${sample}.pe.R2.sai ../bwa_out/${sample}.pe.sam ../bwa_out/${sample}.pe.bam ../bwa_out/${sample}.pe.properly_paired.bam

#        # begin QC
        samtools flagstat ../bwa_out/${sample}.pe.sorted.bam > ../bwa_out/${sample}.flagstat
        mappedReads=\`samtools view  -c ../bwa_out/${sample}.pe.properly_paired.sorted.bam\`
	num_properly_paired=\`grep \"properly paired\" ../bwa_out/${sample}.flagstat | cut -d\" \" -f1\`
	java -jar /share/pkg/picard/1.96/CollectInsertSizeMetrics.jar INPUT=/home/si14w/gnearline/hdc/bwa_out/${sample}.pe.sorted.bam HISTOGRAM_FILE=../results/${sample}.pdf OUTPUT=../results/${sample}.txt
	java -jar /share/pkg/picard/1.96/MarkDuplicates.jar INPUT=/home/si14w/gnearline/hdc/bwa_out/${sample}.pe.properly_paired.sorted.bam METRICS_FILE=../results/${sample}.pcrdups.properly_paired.txt REMOVE_DUPLICATES=true OUTPUT=../bwa_out/${sample}.retest2.deduped.properly_paired.bam
	rm /home/si14w/gnearline/hdc/bwa_out/${sample}.pe.sorted.bam
	percent_dups=\`grep \"Unknown Library\" ../results/${sample}.pcrdups.properly_paired.txt | awk '{ print \$(NF-1)*100}'\`
	totalreads=\$((\${r1_lines} / 4))
	dedupedreads=\`samtools view -c ../bwa_out/${sample}.retest2.deduped.properly_paired.bam\`
        deDupedMappedReadsInPromoter=\`samtools view -F 4 -b ../bwa_out/${sample}.retest2.deduped.properly_paired.bam | bedtools bamtobed -i stdin | bedtools intersect -a stdin -b ../txt/hg19_promoter.bed -u | wc -l\`
	echo ${sample},\${totalreads}, \${mappedReads}, \${dedupedreads},\${deDupedMappedReadsInPromoter},\${percent_dups} > ../txt/${sample}.retest2.stats
        scalingFactor=\`echo \$dedupedreads | awk -vmappedReads=\$dedupedreads 'BEGIN{print 1000000/mappedReads}'\`
        bedtools genomecov -ibam ../bwa_out/${sample}.retest2.deduped.properly_paired.bam -g ~/gnearline/common_scripts/hg19.genome -bga -scale \${scalingFactor} > ../bwa_out/${sample}.bedGraph
        bedGraphToBigWig ../bwa_out/${sample}.bedGraph ~/gnearline/common_scripts/hg19.sorted.genome ../bwa_out/${sample}.bw
	rm ../bwa_out/${sample}.bedGraph
        for exp in {\"high\",\"med\",\"low\",\"verylow\"}
        do
                bedtools slop -i ../txt/gene_tss_expression_time_${time}.\${exp}.minus.bed -b 1000 -g ~/gnearline/common_scripts/hg19.genome | bedtools makewindows -b stdin -w 10 -reverse -i srcwinnum | bigWigAverageOverBed ../bwa_out/${sample}.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/${sample}.\${exp}.minus.bed
                bedtools slop -i ../txt/gene_tss_expression_time_${time}.\${exp}.plus.bed -b 1000 -g ~/gnearline/common_scripts/hg19.genome | bedtools makewindows -b stdin -w 10 -i srcwinnum | bigWigAverageOverBed ../bwa_out/${sample}.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/${sample}.\${exp}.plus.bed
                sed 's/_/\t/g' ../txt/${sample}.\${exp}.minus.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' > ../txt/${sample}.\${exp}.all_signals.txt
                sed 's/_/\t/g' ../txt/${sample}.\${exp}.plus.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' >> ../txt/${sample}.\${exp}.all_signals.txt
        done
        Rscript aggrPlots.R ../txt/${sample}.high.all_signals.txt ../txt/${sample}.med.all_signals.txt ../txt/${sample}.low.all_signals.txt ../txt/${sample}.verylow.all_signals.txt ../results/${sample}.retest2.aggr.pdf $sample
	igvtools count -z 5 -w 25 -e 250 ../bwa_out/${sample}.retest2.deduped.properly_paired.bam  ../bwa_out/${sample}.retest2.pe.properly_paired.tdf  hg19
	rm ../txt/${sample}.\${exp}.minus.bed ../txt/${sample}.\${exp}.plus.bed 
	""" > ../bsubFiles/chip_retest2_anetta_diagen.${sample}.bsub
done<../txt/hdc_patientsamples_anetta_diagen.txt
