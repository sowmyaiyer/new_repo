#	java -jar /share/pkg/picard/1.96/DownsampleSam.jar INPUT=/home/si14w/gnearline/hdc/bwa_out/H3K4me1_200ng_Ch_1ug.pe.properly_paired.sorted.bam  PROBABILITY=0.7 RANDOM_SEED=null OUTPUT=/home/si14w/gnearline/hdc/bwa_out/H3K4me1_200ng_Ch_1ug.pe.properly_paired.downsampled.bam VALIDATION_STRINGENCY=LENIENT
#	samtools sort /home/si14w/gnearline/hdc/bwa_out/H3K4me1_200ng_Ch_1ug.pe.properly_paired.downsampled.bam -o /home/si14w/gnearline/hdc/bwa_out/H3K4me1_200ng_Ch_1ug.pe.properly_paired.sorted.downsampled.bam -T /home/si14w/gnearline/hdc/bwa_out/H3K4me1_200ng_Ch_1ug.pe.properly_paired.sorted
	echo """
#	java -jar /share/pkg/picard/1.96/DownsampleSam.jar INPUT=/home/si14w/gnearline/hdc/bwa_out/H3K4me1_200ng_Ch_1ug.pe.properly_paired.sorted.bam  PROBABILITY=0.15 RANDOM_SEED=null OUTPUT=/home/si14w/gnearline/hdc/bwa_out/H3K4me1_200ng_Ch_1ug.pe.properly_paired.downsampled.bam VALIDATION_STRINGENCY=LENIENT
       samtools sort /home/si14w/gnearline/hdc/bwa_out/H3K4me1_200ng_Ch_1ug.pe.properly_paired.downsampled.bam -o /home/si14w/gnearline/hdc/bwa_out/H3K4me1_200ng_Ch_1ug.pe.properly_paired.sorted.downsampled.bam -T /home/si14w/gnearline/hdc/bwa_out/H3K4me1_200ng_Ch_1ug.pe.properly_paired.sorted
	mappedReads=\`samtools view  ../bwa_out/H3K4me1_200ng_Ch_1ug.pe.properly_paired.sorted.downsampled.bam | wc -l\`
        mappedReadsInPromoter=\`samtools view -F 4 -b ../bwa_out/H3K4me1_200ng_Ch_1ug.pe.properly_paired.sorted.downsampled.bam | bedtools bamtobed -i stdin | bedtools intersect -a stdin -b ../txt/hg19_promoter.bed -u | wc -l\`
	echo downsampled mapped reads = \$mappedReads > ../bwa_out/H3K4me1_200ng_Ch_1ug.downsampled.flagstat
        echo mapped reads in promoter = \$mappedReadsInPromoter >> ../bwa_out/H3K4me1_200ng_Ch_1ug.downsampled.flagstat	
	scalingFactor=\`echo \$mappedReads | awk -vmappedReads=\$mappedReads 'BEGIN{print 1000000/mappedReads}'\`
        bedtools genomecov -ibam ../bwa_out/H3K4me1_200ng_Ch_1ug.pe.properly_paired.sorted.downsampled.bam -g ~/gnearline/common_scripts/hg19.genome -bga -scale \${scalingFactor} > ../bwa_out/H3K4me1_200ng_Ch_1ug.downsampled.bedGraph
        bedGraphToBigWig ../bwa_out/H3K4me1_200ng_Ch_1ug.downsampled.bedGraph ~/gnearline/common_scripts/hg19.sorted.genome ../bwa_out/H3K4me1_200ng_Ch_1ug.downsampled.bw
        rm ../bwa_out/H3K4me1_200ng_Ch_1ug.downsampled.bedGraph
        for exp in {\"high\",\"med\",\"low\",\"verylow\"}
        do
                bedtools slop -i ../txt/gene_tss_expression_Jurkat.\${exp}.minus.bed -b 1000 -g ~/gnearline/common_scripts/hg19.genome | bedtools makewindows -b stdin -w 10 -reverse -i srcwinnum | bigWigAverageOverBed ../bwa_out/H3K4me1_200ng_Ch_1ug.downsampled.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/H3K4me1_200ng_Ch_1ug.downsampled.\${exp}.minus.bed
                bedtools slop -i ../txt/gene_tss_expression_Jurkat.\${exp}.plus.bed -b 1000 -g ~/gnearline/common_scripts/hg19.genome | bedtools makewindows -b stdin -w 10 -i srcwinnum | bigWigAverageOverBed ../bwa_out/H3K4me1_200ng_Ch_1ug.downsampled.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/H3K4me1_200ng_Ch_1ug.downsampled.\${exp}.plus.bed
                sed 's/_/\t/g' ../txt/H3K4me1_200ng_Ch_1ug.downsampled.\${exp}.minus.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' > ../txt/H3K4me1_200ng_Ch_1ug.downsampled.\${exp}.all_signals.txt
                sed 's/_/\t/g' ../txt/H3K4me1_200ng_Ch_1ug.downsampled.\${exp}.plus.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' >> ../txt/H3K4me1_200ng_Ch_1ug.downsampled.\${exp}.all_signals.txt
        done
        Rscript aggrPlots.R ../txt/H3K4me1_200ng_Ch_1ug.downsampled.high.all_signals.txt ../txt/H3K4me1_200ng_Ch_1ug.downsampled.med.all_signals.txt ../txt/H3K4me1_200ng_Ch_1ug.downsampled.low.all_signals.txt ../txt/H3K4me1_200ng_Ch_1ug.downsampled.verylow.all_signals.txt ../results/H3K4me1_200ng_Ch_1ug.downsampled.aggr.pdf
""" > ../bsubFiles/test_Jurkat_downsample_K4me1.bsub
