	sample="natoli_macrophage_SRR502912"
	echo $sample
	expression_for_celltype="mouse_gene_tss_expression_time_mouse_macrophage_0h"
	echo """
	module load IGVTools/2.3.31
	java -jar /share/pkg/picard/1.96/DownsampleSam.jar INPUT=../bwa_out/${sample}.se.bam  PROBABILITY=0.15 RANDOM_SEED=null OUTPUT=../bwa_out/${sample}.se.downsampled.bam VALIDATION_STRINGENCY=LENIENT
	samtools sort -T ../bwa_out/${sample}.se.downsampled.sorted -o ../bwa_out/${sample}.se.downsampled.sorted.bam ../bwa_out/${sample}.se.downsampled.bam
        samtools index ../bwa_out/${sample}.se.downsampled.sorted.bam

#        # begin QC
        samtools flagstat ../bwa_out/${sample}.se.downsampled.sorted.bam > ../bwa_out/${sample}.downsampled.flagstat
        mappedReads=\`samtools view  -c ../bwa_out/${sample}.se.downsampled.sorted.bam\`
        mappedReadsInPromoter=\`samtools view -F 4 -b ../bwa_out/${sample}.se.downsampled.sorted.bam | bedtools bamtobed -i stdin | bedtools intersect -a stdin -b ../txt/mm9_promoter.bed -u | wc -l\`
        echo mapped reads in promoter = \$mappedReadsInPromoter >> ../bwa_out/${sample}.downsampled.flagstat
	java -jar /share/pkg/picard/1.96/CollectInsertSizeMetrics.jar INPUT=../bwa_out/${sample}.se.downsampled.sorted.bam HISTOGRAM_FILE=../results/${sample}.pdf OUTPUT=../results/${sample}.txt
	totalreads=\$((\${r1_lines} / 4))
	echo ${sample},\${totalreads}, \${mappedReads}, \${mappedReadsInPromoter} > ../txt/${sample}.stats
        scalingFactor=\`echo \$mappedReads | awk -vmappedReads=\$mappedReads 'BEGIN{print 1000000/mappedReads}'\`
        bedtools genomecov -ibam ../bwa_out/${sample}.se.downsampled.sorted.bam -g ~/gnearline/mouse_DC_hot/txt/mm9.genome -bga -scale \${scalingFactor} > ../bwa_out/${sample}.bedGraph
        bedGraphToBigWig ../bwa_out/${sample}.bedGraph ~/gnearline/mouse_DC_hot/txt/mm9.genome ../bwa_out/${sample}.bw
	rm ../bwa_out/${sample}.bedGraph
        for exp in {\"high\",\"med\",\"low\",\"verylow\"}
        do
                bedtools slop -i ../txt/${expression_for_celltype}.\${exp}.minus.bed -b 1000 -g ~/gnearline/mouse_DC_hot/txt/mm9.genome | bedtools makewindows -b stdin -w 10 -reverse -i srcwinnum | bigWigAverageOverBed ../bwa_out/${sample}.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/${sample}.\${exp}.minus.bed
                bedtools slop -i ../txt/${expression_for_celltype}.\${exp}.plus.bed -b 1000 -g ~/gnearline/mouse_DC_hot/txt/mm9.genome | bedtools makewindows -b stdin -w 10 -i srcwinnum | bigWigAverageOverBed ../bwa_out/${sample}.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/${sample}.\${exp}.plus.bed
                sed 's/_/\t/g' ../txt/${sample}.\${exp}.minus.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' > ../txt/${sample}.\${exp}.all_signals.txt
                sed 's/_/\t/g' ../txt/${sample}.\${exp}.plus.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' >> ../txt/${sample}.\${exp}.all_signals.txt
        done
        Rscript aggrPlots.R ../txt/${sample}.high.all_signals.txt ../txt/${sample}.med.all_signals.txt ../txt/${sample}.low.all_signals.txt ../txt/${sample}.verylow.all_signals.txt ../results/${sample}.aggr.pdf $sample
	igvtools count -z 5 -w 25 -e 200 ../bwa_out/${sample}.se.sorted.bam  ../bwa_out/${sample}.se.sorted.tdf  mm9
	""" > ../bsubFiles/mouse_macrophage_natoli_downsampled_${sample}.bsub
