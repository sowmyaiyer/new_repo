	echo """
	java -jar /share/pkg/picard/1.96/DownsampleSam.jar INPUT=~/gnearline/hdc/bwa_out/mouse_2h.bam  PROBABILITY=0.15 RANDOM_SEED=null OUTPUT=~/gnearline/hdc/bwa_out/mouse_2h.downsampled.bam VALIDATION_STRINGENCY=LENIENT
	java -jar /share/pkg/picard/1.96/DownsampleSam.jar INPUT=~/gnearline/hdc/bwa_out/mouse_0h.bam  PROBABILITY=0.15 RANDOM_SEED=null OUTPUT=~/gnearline/hdc/bwa_out/mouse_0h.downsampled.bam VALIDATION_STRINGENCY=LENIENT
	samtools merge ~/gnearline/hdc/bwa_out/mouse_0h_and_2h.downsampled.merged.bam ~/gnearline/hdc/bwa_out/mouse_0h.downsampled.bam ~/gnearline/hdc/bwa_out/mouse_2h.downsampled.bam
       samtools sort ~/gnearline/hdc/bwa_out/mouse_0h_and_2h.downsampled.merged.bam -o ~/gnearline/hdc/bwa_out/mouse_0h_and_2h.downsampled.merged.sorted.bam -T ~/gnearline/hdc/bwa_out/mouse_0h_and_2h.downsampled.merged.sorted
	mappedReads=\`samtools view  ~/gnearline/hdc/bwa_out/mouse_0h_and_2h.downsampled.merged.sorted.bam | wc -l\`
        mappedReadsInPromoter=\`samtools view -b ~/gnearline/hdc/bwa_out/mouse_0h_and_2h.downsampled.merged.sorted.bam | bedtools bamtobed -i stdin | bedtools intersect -a stdin -b ../txt/mm9_promoter.bed -u | wc -l\`
	scalingFactor=\`echo \$mappedReads | awk -vmappedReads=\$mappedReads 'BEGIN{print 1000000/mappedReads}'\`
        bedtools genomecov -ibam ~/gnearline/hdc/bwa_out/mouse_0h_and_2h.downsampled.merged.sorted.bam -g /home/si14w/gnearline/mouse_DC_hot/txt/mm9.genome -bga -scale \${scalingFactor} > ../bwa_out/mouse_0h_and_2h.downsampled.merged.sorted.bedGraph
        bedGraphToBigWig ../bwa_out/mouse_0h_and_2h.downsampled.merged.sorted.bedGraph /home/si14w/gnearline/mouse_DC_hot/txt/mm9.genome ../bwa_out/mouse_0h_and_2h.downsampled.merged.sorted.bw
        rm ../bwa_out/mouse_0h_and_2h.downsampled.merged.sorted.bedGraph
        echo \"mouse\",\${mappedReads}, \${mappedReadsInPromoter},\${percent_dups} > ../txt/mouse_H3K27ac_0h.stats
        for exp in {\"high\",\"med\",\"low\",\"verylow\"}
        do
                bedtools slop -i ../txt/mouse_gene_tss_expression_time_mouse_0h.\${exp}.minus.bed -b 1000 -g /home/si14w/gnearline/mouse_DC_hot/txt/mm9.genome | bedtools makewindows -b stdin -w 10 -reverse -i srcwinnum | bigWigAverageOverBed ../bwa_out/mouse_0h_and_2h.downsampled.merged.sorted.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/mouse_0h_and_2h.downsampled.merged.sorted.minus.bed
                bedtools slop -i ../txt/mouse_gene_tss_expression_time_mouse_0h.\${exp}.plus.bed -b 1000 -g /home/si14w/gnearline/mouse_DC_hot/txt/mm9.genome | bedtools makewindows -b stdin -w 10 -i srcwinnum | bigWigAverageOverBed ../bwa_out/mouse_0h_and_2h.downsampled.merged.sorted.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/mouse_0h_and_2h.downsampled.merged.sorted.plus.bed
                sed 's/_/\t/g' ../txt/mouse_0h_and_2h.downsampled.merged.sorted.minus.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' > ../txt/mouse_0h_and_2h.downsampled.merged.\${exp}.all_signals.txt
                sed 's/_/\t/g' ../txt/mouse_0h_and_2h.downsampled.merged.sorted.plus.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' >>  ../txt/mouse_0h_and_2h.downsampled.merged.\${exp}.all_signals.txt
        done
        Rscript aggrPlots.R ../txt/mouse_0h_and_2h.downsampled.merged.high.all_signals.txt ../txt/mouse_0h_and_2h.downsampled.merged.med.all_signals.txt ../txt/mouse_0h_and_2h.downsampled.merged.low.all_signals.txt ../txt/mouse_0h_and_2h.downsampled.merged.verylow.all_signals.txt ../results/mouse_H3K27Ac_0h.aggr.pdf mouse_H3K27A_downsampled_merged
""" > ../bsubFiles/test_mouse_DC_H3K27ac_downsampled.bsub
