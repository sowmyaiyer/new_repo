# get all TSSs
awk '{ if ($6 == "+") tss=$7; else if ($6 == "-") tss=$8; print $1"\t"tss"\t"tss"\t"$4"\t.\t"$6 }' /project/umw_garberlab/edonnard/dc_project/intron_analysis/mouse/reference/mm9_one_tr_per_gene.bed | sort -k1,1 -k2,2n > ../txt/mouse_gene_tss.bed
sort -k4,4 ../txt/mouse_gene_tss.bed > ../txt/mouse_gene_tss.sorted_by_name.bed
# got TPMS from Pranitha (processed RNA seq from Natoli macrophage paper)
awk '{ print $1"\t"$6}' /farline/umw_garberlab/mouse/macropages/rsem_natoli/rsem_macs.genes.results | sort -k1,1 > ~/gnearline/hdc/txt/macrophages.tpm.tsv

# join TSS with gene expression tables
	time="mouse_macrophage_0h"
	join ../txt/mouse_gene_tss.sorted_by_name.bed ~/gnearline/hdc/txt/macrophages.tpm.tsv -1 4 -2 1 | awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$3,$4,$1,$5,$6,$7}' | sort -k1,1 -k2,2n > ../txt/mouse_gene_tss_expression_time_${time}.bed
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
	
