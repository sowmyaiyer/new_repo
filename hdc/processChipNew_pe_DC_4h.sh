for sample in  `awk '{ print $1}' ../txt/chip1222_sampleInfo_uniq.txt`
do
	echo $sample
        for exp in {"high","med","low","verylow"}
        do
                bedtools slop -i ../txt/gene_tss_expression_time_4h.${exp}.minus.bed -b 1000 -g ~/gnearline/common_scripts/hg19.genome | bedtools makewindows -b stdin -w 10 -reverse -i srcwinnum | bigWigAverageOverBed ../bwa_out/${sample}.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/${sample}.${exp}.minus.bed
                bedtools slop -i ../txt/gene_tss_expression_time_4h.${exp}.plus.bed -b 1000 -g ~/gnearline/common_scripts/hg19.genome | bedtools makewindows -b stdin -w 10 -i srcwinnum | bigWigAverageOverBed ../bwa_out/${sample}.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/${sample}.${exp}.plus.bed
                sed 's/_/\t/g' ../txt/${sample}.${exp}.minus.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' > ../txt/${sample}.${exp}.all_signals.txt
                sed 's/_/\t/g' ../txt/${sample}.${exp}.plus.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' >> ../txt/${sample}.${exp}.all_signals.txt
        done
        Rscript aggrPlots.R ../txt/${sample}.high.all_signals.txt ../txt/${sample}.med.all_signals.txt ../txt/${sample}.low.all_signals.txt ../txt/${sample}.verylow.all_signals.txt ../results/${sample}.aggr.pdf
	igvtools count -z 5 -w 25 -e 250 ../bwa_out/${sample}.pe.properly_paired.sorted.bam  ../bwa_out/${sample}.pe.properly_paired.tdf  hg19
done
