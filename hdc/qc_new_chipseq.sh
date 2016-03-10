for sample in  `awk '{ print $1}' ../txt/chipnew_sampleInfo.txt | sort | uniq`
do
        echo $sample

        #samtools flagstat ../bwa_out/${sample}.se.sorted.bam > ../bwa_out/${sample}.se.flagstat
        #samtools flagstat ../bwa_out/${sample}.pe.sorted.bam > ../bwa_out/${sample}.pe.flagstat
        
	mappedReads_se=`awk '{ if (NR ==5) print $1}' ../bwa_out/${sample}.se.flagstat`
        mappedReads_pe=`awk '{ if (NR ==5) print $1}' ../bwa_out/${sample}.pe.flagstat`
        
	#mappedReadsInPromoter_pe=`samtools view -F 0x4 -b ../bwa_out/${sample}.pe.sorted.bam | bedtools bamtobed -i stdin | bedtools intersect -a stdin -b ../txt/hg19_promoter.bed -u | wc -l`
        #mappedReadsInPromoter_se=`samtools view -F 0x4 -b ../bwa_out/${sample}.se.sorted.bam | bedtools bamtobed -i stdin | bedtools intersect -a stdin -b ../txt/hg19_promoter.bed -u | wc -l`
        
	#echo mapped reads in promoter = $mappedReadsInPromoter_pe >> ../bwa_out/${sample}.pe.flagstat
        #echo mapped reads in promoter = $mappedReadsInPromoter_se >> ../bwa_out/${sample}.se.flagstat


        scalingFactor_se=`echo $mappedReads_se | awk -vmappedReads=$mappedReads_se 'BEGIN{print 1000000/mappedReads}'`
        scalingFactor_pe=`echo $mappedReads_pe | awk -vmappedReads=$mappedReads_pe 'BEGIN{print 1000000/mappedReads}'`

	echo $scalingFactor_se $scalingFactor_pe

        bedtools genomecov -ibam ../bwa_out/${sample}.se.sorted.bam -g ~/gnearline/common_scripts/hg19.genome -bga -scale ${scalingFactor_se} > ../bwa_out/${sample}.se.bedGraph
        bedtools genomecov -ibam ../bwa_out/${sample}.pe.sorted.bam -g ~/gnearline/common_scripts/hg19.genome -bga -scale ${scalingFactor_pe} > ../bwa_out/${sample}.pe.bedGraph

        bedGraphToBigWig ../bwa_out/${sample}.se.bedGraph ~/gnearline/common_scripts/hg19.genome ../bwa_out/${sample}.se.bw
        bedGraphToBigWig ../bwa_out/${sample}.pe.bedGraph ~/gnearline/common_scripts/hg19.genome ../bwa_out/${sample}.pe.bw
	
	time="0h"
        for exp in {"high","med","low","verylow"}
        do
                bedtools slop -i ../txt/gene_tss_expression_time_${time}.${exp}.minus.bed -b 1000 -g ~/gnearline/common_scripts/hg19.genome | bedtools makewindows -b stdin -w 10 -reverse -i srcwinnum | bigWigAverageOverBed ../bwa_out/${sample}.se.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/${sample}.${exp}.se.minus.bed
                bedtools slop -i ../txt/gene_tss_expression_time_${time}.${exp}.plus.bed -b 1000 -g ~/gnearline/common_scripts/hg19.genome | bedtools makewindows -b stdin -w 10 -i srcwinnum | bigWigAverageOverBed ../bwa_out/${sample}.se.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/${sample}.${exp}.se.plus.bed
                bedtools slop -i ../txt/gene_tss_expression_time_${time}.${exp}.minus.bed -b 1000 -g ~/gnearline/common_scripts/hg19.genome | bedtools makewindows -b stdin -w 10 -reverse -i srcwinnum | bigWigAverageOverBed ../bwa_out/${sample}.pe.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/${sample}.${exp}.pe.minus.bed
                bedtools slop -i ../txt/gene_tss_expression_time_${time}.${exp}.plus.bed -b 1000 -g ~/gnearline/common_scripts/hg19.genome | bedtools makewindows -b stdin -w 10 -i srcwinnum | bigWigAverageOverBed ../bwa_out/${sample}.pe.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/${sample}.${exp}.pe.plus.bed
                sed 's/_/\t/g' ../txt/${sample}.${exp}.se.minus.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' > ../txt/${sample}.${exp}.se.all_signals.txt
                sed 's/_/\t/g' ../txt/${sample}.${exp}.se.plus.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' >> ../txt/${sample}.${exp}.se.all_signals.txt
                sed 's/_/\t/g' ../txt/${sample}.${exp}.pe.minus.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' > ../txt/${sample}.${exp}.pe.all_signals.txt
                sed 's/_/\t/g' ../txt/${sample}.${exp}.pe.plus.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' >> ../txt/${sample}.${exp}.pe.all_signals.txt
        done
        sed -i '/GTF2H2C/d' ../txt/${sample}.med.se.all_signals.txt
        sed -i '/GTF2H2C/d' ../txt/${sample}.med.pe.all_signals.txt
        Rscript aggrPlots.R ../txt/${sample}.high.se.all_signals.txt ../txt/${sample}.med.se.all_signals.txt ../txt/${sample}.low.se.all_signals.txt ../txt/${sample}.verylow.se.all_signals.txt ../results/${sample}.H3K4me3.se.pdf ${sample}
        Rscript aggrPlots.R ../txt/${sample}.high.pe.all_signals.txt ../txt/${sample}.med.pe.all_signals.txt ../txt/${sample}.low.pe.all_signals.txt ../txt/${sample}.verylow.pe.all_signals.txt ../results/${sample}.H3K4me3.pe.pdf ${sample}
done
