# ./processChipNew_pe_JK_vs_DC_preprocess.sh
while read sample
do
	time_tmp=`echo $sample | cut -d"_" -f2`
	time_tmp2=${time_tmp}"_LPS"
	if [[ ($time_tmp == "CV") || ($time_tmp == "PM") || ($time_tmp == "30m") ]]; then
		time="0h"
	else
		time=$time_tmp2
	fi
	echo $time_tmp2 $time
	echo """
	for exp in {\"high\",\"med\",\"low\",\"verylow\"}
        do
		echo start \${exp}
                bedtools slop -i ../txt/gene_tss_expression_time_${time}.\${exp}.minus.bed -b 1000 -g ~/gnearline/common_scripts/hg19.genome | bedtools makewindows -b stdin -w 10 -reverse -i srcwinnum | bigWigAverageOverBed ../bwa_out/${sample}.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/${sample}.\${exp}.minus.bed
                bedtools slop -i ../txt/gene_tss_expression_time_${time}.\${exp}.plus.bed -b 1000 -g ~/gnearline/common_scripts/hg19.genome | bedtools makewindows -b stdin -w 10 -i srcwinnum | bigWigAverageOverBed ../bwa_out/${sample}.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/${sample}.\${exp}.plus.bed
                sed 's/_/\t/g' ../txt/${sample}.\${exp}.minus.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' > ../txt/${sample}.\${exp}.all_signals.txt
                sed 's/_/\t/g' ../txt/${sample}.\${exp}.plus.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' >> ../txt/${sample}.\${exp}.all_signals.txt
		echo done \${exp}
        done
        Rscript aggrPlots.R $HOME/gnearline/hdc/txt/${sample}.high.all_signals.txt $HOME/gnearline/hdc/txt/${sample}.med.all_signals.txt $HOME/gnearline/hdc/txt/${sample}.low.all_signals.txt $HOME/gnearline/hdc/txt/${sample}.verylow.all_signals.txt $HOME/gnearline/hdc/results/${sample}.aggr.pdf $sample """ > ../bsubFiles/aggrplots_newrnaseq.${sample}.bsub
done<../txt/hdc_patientsamples_info_0205.sub.txt
