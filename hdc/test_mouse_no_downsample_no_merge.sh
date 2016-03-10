	echo """
	mappedReads=\`samtools view  ~/gnearline/hdc/bwa_out/mouse_0h_K4me1.bam | wc -l\`
        mappedReadsInPromoter=\`samtools view -b ~/gnearline/hdc/bwa_out/mouse_0h_K4me1.bam | bedtools bamtobed -i stdin | bedtools intersect -a stdin -b ../txt/mm9_promoter.bed -u | wc -l\`
	scalingFactor=\`echo \$mappedReads | awk -vmappedReads=\$mappedReads 'BEGIN{print 1000000/mappedReads}'\`
        bedtools genomecov -ibam ~/gnearline/hdc/bwa_out/mouse_0h_K4me1.bam -g /home/si14w/gnearline/mouse_DC_hot/txt/mm9.genome -bga -scale \${scalingFactor} > ../bwa_out/mouse_0h.K4me1.bedGraph
        bedGraphToBigWig ../bwa_out/mouse_0h.K4me1.bedGraph /home/si14w/gnearline/mouse_DC_hot/txt/mm9.genome ../bwa_out/mouse_0h.K4me1.bw
        rm ../bwa_out/mouse_0h.K4me1.bedGraph
        echo \"mouse\",\${mappedReads}, \${mappedReadsInPromoter} > ../txt/mouse_H3K4me1_0h_nodownsample.stats
        for exp in {\"high\",\"med\",\"low\",\"verylow\"}
        do
                bedtools slop -i ../txt/mouse_gene_tss_expression_time_mouse_0h.\${exp}.minus.bed -b 1000 -g /home/si14w/gnearline/mouse_DC_hot/txt/mm9.genome | bedtools makewindows -b stdin -w 10 -reverse -i srcwinnum | bigWigAverageOverBed ../bwa_out/mouse_0h.K4me1.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/mouse_0h.K4me1.minus.bed
                bedtools slop -i ../txt/mouse_gene_tss_expression_time_mouse_0h.\${exp}.plus.bed -b 1000 -g /home/si14w/gnearline/mouse_DC_hot/txt/mm9.genome | bedtools makewindows -b stdin -w 10 -i srcwinnum | bigWigAverageOverBed ../bwa_out/mouse_0h.K4me1.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/mouse_0h.K4me1.plus.bed
                sed 's/_/\t/g' ../txt/mouse_0h.K4me1.minus.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' > ../txt/mouse_0h.K4me1.\${exp}.all_signals.txt
                sed 's/_/\t/g' ../txt/mouse_0h.K4me1.plus.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' >>  ../txt/mouse_0h.K4me1.\${exp}.all_signals.txt
        done
        Rscript aggrPlots.R ../txt/mouse_0h.K4me1.high.all_signals.txt ../txt/mouse_0h.K4me1.med.all_signals.txt ../txt/mouse_0h.K4me1.low.all_signals.txt ../txt/mouse_0h.K4me1.verylow.all_signals.txt ../results/mouse_H3K4me1_0h.aggr.pdf mouse_H3K4me1_0h_nodownsample
""" > ../bsubFiles/test_mouse_DC_H3K4me1_nodownsampled.bsub
