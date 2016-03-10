# get all TSSs
# cp /share/data/umw_biocore/genome_data/rsem_correction/refseq_hg19/hg19_one_tr_per_gene.bed ../txt/gene_tss.bed
awk '{ if ($6 == "+") tss=$7; else if ($6 == "-") tss=$8; print $1"\t"tss"\t"tss"\t"$4"\t.\t"$6 }' /project/umw_garberlab/edonnard/dc_project/intron_analysis/mouse/reference/mm9_one_tr_per_gene.bed | sort -k1,1 -k2,2n > ../txt/mouse_gene_tss.bed
sort -k4,4 ../txt/mouse_gene_tss.bed > ../txt/mouse_gene_tss.sorted_by_name.bed
# get all expression values
awk '{ if (NR > 1) print $1"\t"$3}' ~/gnearline/mouse_DC_hot/txt/Mouse_tpm.tsv > ../txt/mouse_0h_tpms.txt

# join TSS with gene expression tables
	time="mouse_0h"
	join ../txt/mouse_gene_tss.sorted_by_name.bed ../txt/mouse_0h_tpms.txt -1 4 -2 1 | awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$3,$4,$1,$5,$6,$7}' | sort -k1,1 -k2,2n > ../txt/mouse_gene_tss_expression_time_${time}.bed
	quartiles=`awk '{ print $NF}' ../txt/mouse_gene_tss_expression_time_${time}.bed | Rscript getQuantile.R`
	q1=`echo $quartiles | cut -d" " -f4`
	q2=`echo $quartiles | cut -d" " -f3`
	q3=`echo $quartiles | cut -d" " -f2`
	awk -vq1=$q1 '{ if ($NF > q1 && $6 == "-") print }' ../txt/mouse_gene_tss_expression_time_${time}.bed > ../txt/mouse_gene_tss_expression_time_${time}.high.minus.bed
	awk -vq1=$q1 '{ if ($NF > q1 && $6 == "+") print }' ../txt/mouse_gene_tss_expression_time_${time}.bed > ../txt/mouse_gene_tss_expression_time_${time}.high.plus.bed

	awk -vq1=$q1 -vq2=$q2 '{ if ($NF <= q1 && $NF > q2 && $6 == "-") print }' ../txt/mouse_gene_tss_expression_time_${time}.bed > ../txt/mouse_gene_tss_expression_time_${time}.med.minus.bed
	awk -vq1=$q1 -vq2=$q2 '{ if ($NF <= q1 && $NF > q2 && $6 == "+") print }' ../txt/mouse_gene_tss_expression_time_${time}.bed > ../txt/mouse_gene_tss_expression_time_${time}.med.plus.bed

	awk -vq2=$q2 -vq3=$q3 '{ if ($NF <= q2 && $NF > q3 && $6 == "-") print }' ../txt/mouse_gene_tss_expression_time_${time}.bed > ../txt/mouse_gene_tss_expression_time_${time}.low.minus.bed
	awk -vq2=$q2 -vq3=$q3 '{ if ($NF <= q2 && $NF > q3 && $6 == "+") print }' ../txt/mouse_gene_tss_expression_time_${time}.bed > ../txt/mouse_gene_tss_expression_time_${time}.low.plus.bed

	awk -vq3=$q3 '{ if ($NF <= q3 && $6 == "-") print }' ../txt/mouse_gene_tss_expression_time_${time}.bed > ../txt/mouse_gene_tss_expression_time_${time}.verylow.minus.bed
	awk -vq3=$q3 '{ if ($NF <= q3 && $6 == "+") print }' ../txt/mouse_gene_tss_expression_time_${time}.bed > ../txt/mouse_gene_tss_expression_time_${time}.verylow.plus.bed
	
