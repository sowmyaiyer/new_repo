# Got ensemble gene coords from UCSC table browser
#awk 'BEGIN{OFS="\t"}{print $2,$4,$5,$1,".",$3}' ../txt/ensembl.gene-3.bed | sort -k4,4 | bedtools groupby -g 4 -c 1,2,3 -o first,first,first -full | awk 'BEGIN{OFS="\t"}{ if ($6 == "+") tss=$2; else if ($6 == "-") tss=$3; print $1"\t"tss"\t"tss"\t"$4"\t.\t"$6 }' | sort -k4,4 > ../txt/ensembl.tss.sorted_by_transcriptid.bed
#zcat ../chipJurkat/ENCFF000MSP.gtf.gz | awk -F"\t" '{ if ($3 == "transcript") {split($9,arr," ");print arr[4]"\t"arr[6]}}'  | sed 's/\"//g' | sed 's/;//g' | sort -k1,1 | uniq > ../txt/Jurkat_rnaseq_fpkm.txt

#join ../txt/ensembl.tss.sorted_by_transcriptid.bed ../txt/Jurkat_rnaseq_fpkm.txt -1 4 -2 1 | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$1,$5,$6,$7}' |  sort -k1,1 -k2,2n > ../txt/gene_tss_expression_Jurkat.bed
quartiles=`awk '{ print $NF}' ../txt/gene_tss_expression_Jurkat.bed | Rscript getQuantile.R`
q1=`echo $quartiles | cut -d" " -f4`
q2=`echo $quartiles | cut -d" " -f3`
q3=`echo $quartiles | cut -d" " -f2`
awk -vq1=$q1 '{ if ($NF > q1 && $6 == "-") print }' ../txt/gene_tss_expression_Jurkat.bed > ../txt/gene_tss_expression_Jurkat.high.minus.bed
awk -vq1=$q1 '{ if ($NF > q1 && $6 == "+") print }' ../txt/gene_tss_expression_Jurkat.bed > ../txt/gene_tss_expression_Jurkat.high.plus.bed

awk -vq1=$q1 -vq2=$q2 '{ if ($NF <= q1 && $NF > q2 && $6 == "-") print }' ../txt/gene_tss_expression_Jurkat.bed > ../txt/gene_tss_expression_Jurkat.med.minus.bed
awk -vq1=$q1 -vq2=$q2 '{ if ($NF <= q1 && $NF > q2 && $6 == "+") print }' ../txt/gene_tss_expression_Jurkat.bed > ../txt/gene_tss_expression_Jurkat.med.plus.bed

awk -vq2=$q2 -vq3=$q3 '{ if ($NF <= q2 && $NF > q3 && $6 == "-") print }' ../txt/gene_tss_expression_Jurkat.bed > ../txt/gene_tss_expression_Jurkat.low.minus.bed
awk -vq2=$q2 -vq3=$q3 '{ if ($NF <= q2 && $NF > q3 && $6 == "+") print }' ../txt/gene_tss_expression_Jurkat.bed > ../txt/gene_tss_expression_Jurkat.low.plus.bed

awk -vq3=$q3 '{ if ($NF <= q3 && $6 == "-") print }' ../txt/gene_tss_expression_Jurkat.bed > ../txt/gene_tss_expression_Jurkat.verylow.minus.bed
awk -vq3=$q3 '{ if ($NF <= q3 && $6 == "+") print }' ../txt/gene_tss_expression_Jurkat.bed > ../txt/gene_tss_expression_Jurkat.verylow.plus.bed

for exp in {"high","med","low","verylow"}
do
		echo $exp
        	bedtools slop -i ../txt/gene_tss_expression_Jurkat.${exp}.minus.bed -b 1000 -g ~/gnearline/common_scripts/hg19.genome | bedtools makewindows -b stdin -w 10 -reverse -i srcwinnum | bigWigAverageOverBed ../chipJurkat/ENCFF001FTP_Jurkat_H3K4me3_ENCODE.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/ENCFF001FTP.Jurkat.H3K4me3.minus.${exp}.bed
                bedtools slop -i ../txt/gene_tss_expression_Jurkat.${exp}.plus.bed -b 1000 -g ~/gnearline/common_scripts/hg19.genome | bedtools makewindows -b stdin -w 10 -i srcwinnum | bigWigAverageOverBed ../chipJurkat/ENCFF001FTP_Jurkat_H3K4me3_ENCODE.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/ENCFF001FTP.Jurkat.H3K4me3.plus.${exp}.bed
                sed 's/_/\t/g' ../txt/ENCFF001FTP.Jurkat.H3K4me3.minus.${exp}.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' > ../txt/ENCFF001FTP.Jurkat.H3K4me3.${exp}.all_signals.txt
                sed 's/_/\t/g' ../txt/ENCFF001FTP.Jurkat.H3K4me3.plus.${exp}.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' >> ../txt/ENCFF001FTP.Jurkat.H3K4me3.${exp}.all_signals.txt
done
