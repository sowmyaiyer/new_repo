# get all TSSs
# cp /share/data/umw_biocore/genome_data/rsem_correction/refseq_hg19/hg19_one_tr_per_gene.bed ../txt/gene_tss.bed
awk '{ if ($6 == "+") tss=$7; else if ($6 == "-") tss=$8; print $1"\t"tss"\t"tss"\t"$4"\t.\t"$6 }' ../txt/gene_tss.bed | sort -k1,1 -k2,2n > ../txt/gene_tss_sorted.bed
sort -k4,4 ../txt/gene_tss_sorted.bed > ../txt/gene_tss_sorted_by_name.bed
# get all expression values
awk '{ if(NR > 1) print $1"\t"$2}' ../txt/D50_6d_LPS.tsv | sort -k1,1 > ../txt/D50_6d_LPS_time_0h.tsv
awk '{ if(NR > 1) print $1"\t"$3}' ../txt/D50_6d_LPS.tsv | sort -k1,1 > ../txt/D50_6d_LPS_time_2h_LPS.tsv
awk '{ if(NR > 1) print $1"\t"$4}' ../txt/D50_6d_LPS.tsv | sort -k1,1 > ../txt/D50_6d_LPS_time_4h_LPS.tsv

# join TSS with gene expression tables
for time in {"0h","2h_LPS","4h_LPS"}
do
	join ../txt/gene_tss_sorted_by_name.bed ../txt/D50_6d_LPS_time_${time}.tsv -1 4 -2 1 | awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$3,$4,$1,$5,$6,$7}' | sort -k1,1 -k2,2n > ../txt/gene_tss_expression_time_${time}.bed
	quartiles=`awk '{ print $NF}' ../txt/gene_tss_expression_time_${time}.bed | Rscript getQuantile.R`
	q1=`echo $quartiles | cut -d" " -f4`
	q2=`echo $quartiles | cut -d" " -f3`
	q3=`echo $quartiles | cut -d" " -f2`
	awk -vq1=$q1 '{ if ($NF > q1 && $6 == "-") print }' ../txt/gene_tss_expression_time_${time}.bed > ../txt/gene_tss_expression_time_${time}.high.minus.bed
	awk -vq1=$q1 '{ if ($NF > q1 && $6 == "+") print }' ../txt/gene_tss_expression_time_${time}.bed > ../txt/gene_tss_expression_time_${time}.high.plus.bed

	awk -vq1=$q1 -vq2=$q2 '{ if ($NF <= q1 && $NF > q2 && $6 == "-") print }' ../txt/gene_tss_expression_time_${time}.bed > ../txt/gene_tss_expression_time_${time}.med.minus.bed
	awk -vq1=$q1 -vq2=$q2 '{ if ($NF <= q1 && $NF > q2 && $6 == "+") print }' ../txt/gene_tss_expression_time_${time}.bed > ../txt/gene_tss_expression_time_${time}.med.plus.bed

	awk -vq2=$q2 -vq3=$q3 '{ if ($NF <= q2 && $NF > q3 && $6 == "-") print }' ../txt/gene_tss_expression_time_${time}.bed > ../txt/gene_tss_expression_time_${time}.low.minus.bed
	awk -vq2=$q2 -vq3=$q3 '{ if ($NF <= q2 && $NF > q3 && $6 == "+") print }' ../txt/gene_tss_expression_time_${time}.bed > ../txt/gene_tss_expression_time_${time}.low.plus.bed

	awk -vq3=$q3 '{ if ($NF <= q3 && $6 == "-") print }' ../txt/gene_tss_expression_time_${time}.bed > ../txt/gene_tss_expression_time_${time}.verylow.minus.bed
	awk -vq3=$q3 '{ if ($NF <= q3 && $6 == "+") print }' ../txt/gene_tss_expression_time_${time}.bed > ../txt/gene_tss_expression_time_${time}.verylow.plus.bed
	
	for exp in {"high","med","low","verylow"}
	do
		for ab in {"H3K4me3","H3K27ac","H3K4me1"}
		do
		bedtools slop -i ../txt/gene_tss_expression_time_${time}.${exp}.minus.bed -b 1000 -g ~/gnearline/common_scripts/hg19.genome | bedtools makewindows -b stdin -w 10 -reverse -i srcwinnum | bigWigAverageOverBed ../chip/ucsc_chip/E01_6d_${time}_${ab}.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/E01_6d_${time}_${ab}.${time}.${exp}.minus.bed
		bedtools slop -i ../txt/gene_tss_expression_time_${time}.${exp}.plus.bed -b 1000 -g ~/gnearline/common_scripts/hg19.genome | bedtools makewindows -b stdin -w 10 -i srcwinnum | bigWigAverageOverBed ../chip/ucsc_chip/E01_6d_${time}_${ab}.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/E01_6d_${time}_${ab}.${time}.${exp}.plus.bed
		sed 's/_/\t/g' ../txt/E01_6d_${time}_${ab}.${time}.${exp}.minus.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' > ../txt/E01_6d_${time}_${ab}.${time}.${exp}.all_signals.txt	
		sed 's/_/\t/g' ../txt/E01_6d_${time}_${ab}.${time}.${exp}.plus.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' >> ../txt/E01_6d_${time}_${ab}.${time}.${exp}.all_signals.txt	
		done
	done
	Rscript aggrPlots.R ../txt/E01_6d_${time}_H3K4me3.${time}.high.all_signals.txt ../txt/E01_6d_${time}_H3K4me3.${time}.med.all_signals.txt ../txt/E01_6d_${time}_H3K4me3.${time}.low.all_signals.txt ../txt/E01_6d_${time}_H3K4me3.${time}.verylow.all_signals.txt ../results/E01_6d_${time}_H3K4me3.aggr.pdf
	Rscript aggrPlots.R ../txt/E01_6d_${time}_H3K27ac.${time}.high.all_signals.txt ../txt/E01_6d_${time}_H3K27ac.${time}.med.all_signals.txt ../txt/E01_6d_${time}_H3K27ac.${time}.low.all_signals.txt ../txt/E01_6d_${time}_H3K27ac.${time}.verylow.all_signals.txt ../results/E01_6d_${time}_H3K27ac.aggr.pdf
	Rscript aggrPlots.R ../txt/E01_6d_${time}_H3K4me1.${time}.high.all_signals.txt ../txt/E01_6d_${time}_H3K4me1.${time}.med.all_signals.txt ../txt/E01_6d_${time}_H3K4me1.${time}.low.all_signals.txt ../txt/E01_6d_${time}_H3K4me1.${time}.verylow.all_signals.txt ../results/E01_6d_${time}_H3K4me1.aggr.pdf
done
