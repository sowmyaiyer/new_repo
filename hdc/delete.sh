while read sample
do
                time="0h"
	echo """
	module load cutadapt/1.9
	module load IGVTools/2.3.31

	dedupedreads=\`samtools view -c ../bwa_out/${sample}.bam\`
        deDupedMappedReadsInPromoter=\`samtools view -F 4 -b ../bwa_out/${sample}.bam | bedtools bamtobed -i stdin | bedtools intersect -a stdin -b ../txt/hg19_promoter.bed -u | wc -l\`
	echo ${sample},\${totalreads}, \${mappedReads}, \${dedupedreads},\${deDupedMappedReadsInPromoter},\${percent_dups} > ../txt/${sample}.retest2.stats
        scalingFactor=\`echo \$dedupedreads | awk -vmappedReads=\$dedupedreads 'BEGIN{print 1000000/mappedReads}'\`
        bedtools genomecov -ibam ../bwa_out/${sample}.bam -g ~/gnearline/common_scripts/hg19.genome -bga -scale \${scalingFactor} > ../bwa_out/${sample}.bedGraph
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
	igvtools count -z 5 -w 25 -e 250 ../bwa_out/${sample}.bam  ../bwa_out/${sample}.retest2.pe.properly_paired.tdf  hg19
	rm ../txt/${sample}.\${exp}.minus.bed ../txt/${sample}.\${exp}.plus.bed 
	""" > ../bsubFiles/chip_retest2_patient_samples_merged.${sample}.bsub
done</home/si14w/gnearline/hdc/txt/mergedSamples.txt
