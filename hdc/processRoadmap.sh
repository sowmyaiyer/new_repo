# Got ensemble gene coords from UCSC table browser
awk 'BEGIN{OFS="\t"}{print $2,$4,$5,$6,".",$3}' ../txt/ensembl.gene-3.bed | sort -k4,4 | bedtools groupby -g 4 -c 1,2,3 -o first,first,first -full | awk 'BEGIN{OFS="\t"}{ if ($6 == "+") tss=$2; else if ($6 == "-") tss=$3; print $1"\t"tss"\t"tss"\t"$4"\t.\t"$6 }' | sort -k4,4 > ../txt/ensembl.tss.sorted_by_geneid.bed
#wget http://egg2.wustl.edu/roadmap/data/byDataType/rna/expression/57epigenomes.RPKM.pc.gz
# Figured out from Anshul's googledoc that E062 is the epigenome ID for this cell type
awk '{print $1"\t"$27}' ../txt/57epigenomes.RPKM.pc | sort -k1,1 > ../txt/RPKM_E062_primary_mononuclear.txt
#wget "http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM1127126&format=file&file=GSM1127126%5FUCSF%2DUBC%2EPeripheral%5FBlood%5FMononuclear%5FPrimary%5FCells%2EH3K4me3%2ETC015%2Ewig%2Egz" -O GSM1127126.roadmap.H3K4me3.gz
#wigToBigWig GSM1127126.roadmap.H3K4me3.wig ~/gnearline/common_scripts/hg19.genome GSM1127126.roadmap.H3K4me3.bw

join ../txt/ensembl.tss.sorted_by_geneid.bed ../txt/RPKM_E062_primary_mononuclear.txt -1 4 -2 1 | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$1,$5,$6,$7}' |  sort -k1,1 -k2,2n > ../txt/gene_tss_expression_roadmap.bed
quartiles=`awk '{ print $NF}' ../txt/gene_tss_expression_roadmap.bed | Rscript getQuantile.R`
q1=`echo $quartiles | cut -d" " -f4`
q2=`echo $quartiles | cut -d" " -f3`
q3=`echo $quartiles | cut -d" " -f2`
awk -vq1=$q1 '{ if ($NF > q1 && $6 == "-") print }' ../txt/gene_tss_expression_roadmap.bed > ../txt/gene_tss_expression_roadmap.high.minus.bed
awk -vq1=$q1 '{ if ($NF > q1 && $6 == "+") print }' ../txt/gene_tss_expression_roadmap.bed > ../txt/gene_tss_expression_roadmap.high.plus.bed

awk -vq1=$q1 -vq2=$q2 '{ if ($NF <= q1 && $NF > q2 && $6 == "-") print }' ../txt/gene_tss_expression_roadmap.bed > ../txt/gene_tss_expression_roadmap.med.minus.bed
awk -vq1=$q1 -vq2=$q2 '{ if ($NF <= q1 && $NF > q2 && $6 == "+") print }' ../txt/gene_tss_expression_roadmap.bed > ../txt/gene_tss_expression_roadmap.med.plus.bed

awk -vq2=$q2 -vq3=$q3 '{ if ($NF <= q2 && $NF > q3 && $6 == "-") print }' ../txt/gene_tss_expression_roadmap.bed > ../txt/gene_tss_expression_roadmap.low.minus.bed
awk -vq2=$q2 -vq3=$q3 '{ if ($NF <= q2 && $NF > q3 && $6 == "+") print }' ../txt/gene_tss_expression_roadmap.bed > ../txt/gene_tss_expression_roadmap.low.plus.bed

awk -vq3=$q3 '{ if ($NF <= q3 && $6 == "-") print }' ../txt/gene_tss_expression_roadmap.bed > ../txt/gene_tss_expression_roadmap.verylow.minus.bed
awk -vq3=$q3 '{ if ($NF <= q3 && $6 == "+") print }' ../txt/gene_tss_expression_roadmap.bed > ../txt/gene_tss_expression_roadmap.verylow.plus.bed

for exp in {"high","med","low","verylow"}
do
		echo $exp
        	bedtools slop -i ../txt/gene_tss_expression_roadmap.${exp}.minus.bed -b 1000 -g ~/gnearline/common_scripts/hg19.genome | bedtools makewindows -b stdin -w 10 -reverse -i srcwinnum | bigWigAverageOverBed ../txt/GSM1127126.roadmap.H3K4me3.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/GSM1127126.roadmap.H3K4me3.minus.${exp}.bed
                bedtools slop -i ../txt/gene_tss_expression_roadmap.${exp}.plus.bed -b 1000 -g ~/gnearline/common_scripts/hg19.genome | bedtools makewindows -b stdin -w 10 -i srcwinnum | bigWigAverageOverBed ../txt/GSM1127126.roadmap.H3K4me3.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/GSM1127126.roadmap.H3K4me3.plus.${exp}.bed
                sed 's/_/\t/g' ../txt/GSM1127126.roadmap.H3K4me3.minus.${exp}.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' > ../txt/GSM1127126.roadmap.H3K4me3.${exp}.all_signals.txt
                sed 's/_/\t/g' ../txt/GSM1127126.roadmap.H3K4me3.plus.${exp}.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' >> ../txt/GSM1127126.roadmap.H3K4me3.${exp}.all_signals.txt
done
