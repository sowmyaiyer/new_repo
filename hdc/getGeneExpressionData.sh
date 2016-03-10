# get all TSSs
# cp /share/data/umw_biocore/genome_data/rsem_correction/refseq_hg19/hg19_one_tr_per_gene.bed ../txt/gene_tss.bed
#awk '{ if ($6 == "+") tss=$7; else if ($6 == "-") tss=$8; print $1"\t"tss"\t"tss"\t"$4"\t.\t"$6 }' ../txt/gene_tss.bed | sort -k1,1 -k2,2n > ../txt/gene_tss_sorted.bed
#sort -k4,4 ../txt/gene_tss_sorted.bed > ../txt/gene_tss_sorted_by_name.bed
# get all expression values
### changed names for 2 genes to change "_" to "." in gene names for two genes APOBEC3A.B and GTF2H2C.2 to accomodate bedtools makewindows function
awk '{ if(NR > 1) print $1"\t"$2}' ../txt/D01_MDDC_LPS.tsv | sort -k1,1 > ../txt/D01_MDDC_LPS_time_0h.tsv
awk '{ if(NR > 1) print $1"\t"$3}' ../txt/D01_MDDC_LPS.tsv | sort -k1,1 > ../txt/D01_MDDC_LPS_time_1h_LPS.tsv
awk '{ if(NR > 1) print $1"\t"$4}' ../txt/D01_MDDC_LPS.tsv | sort -k1,1 > ../txt/D01_MDDC_LPS_time_2h_LPS.tsv
awk '{ if(NR > 1) print $1"\t"$5}' ../txt/D01_MDDC_LPS.tsv | sort -k1,1 > ../txt/D01_MDDC_LPS_time_4h_LPS.tsv
awk '{ if(NR > 1) print $1"\t"$6}' ../txt/D01_MDDC_LPS.tsv | sort -k1,1 > ../txt/D01_MDDC_LPS_time_6h_LPS.tsv
awk '{ if(NR > 1) print $1"\t"$7}' ../txt/D01_MDDC_LPS.tsv | sort -k1,1 > ../txt/D01_MDDC_LPS_time_12h_LPS.tsv
awk '{ if(NR > 1) print $1"\t"$8}' ../txt/D01_MDDC_LPS.tsv | sort -k1,1 > ../txt/D01_MDDC_LPS_time_24h_LPS.tsv

# join TSS with gene expression tables
for time in {"0h","1h_LPS","2h_LPS","4h_LPS","6h_LPS","12h_LPS","24h_LPS"}
do
        join ../txt/gene_tss_sorted_by_name.bed ../txt/D01_MDDC_LPS_time_${time}.tsv -1 4 -2 1 | awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$3,$4,$1,$5,$6,$7}' | sort -k1,1 -k2,2n > ../txt/gene_tss_expression_time_${time}.bed
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
done
