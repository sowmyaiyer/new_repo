## BEGIN expression analysis ###
atac_peak_file=$1
bedtools intersect -a ${atac_peak_file}  -b ../data/H1_and_H7_H3K27ac_merged.bed -u > ../data/H1_ATAC_intersect_H3K27ac.bed
bedtools intersect -a ../data/H1_ATAC.distal.sorted.bed -b ../data/H1_and_H7_H3K27ac_merged.bed -v > ../data/H1_ATAC_non_intersect_H3K27ac.bed
cat ../data/GencodeV19.promoter.bed ../data/H1_ATAC.distal.sorted.bed | awk 'BEGIN{OFS="\t"}{ print $1,$2,$3}' > ../txt/excludable_ATAC_and_promoters.bed

awk '{ if ($3 == "transcript") print }' ../data/H1hesc.transcript.quant.gtf  | awk -f extractFpkmWithTss.awk > ../data/H1hesc.transcript_only.tss_and_fpkm.bed

bedtools closest -a  ../data/H1_ATAC_intersect_H3K27ac.bed -b ../data/H1hesc.transcript_only.tss_and_fpkm.bed | awk 'BEGIN{OFS="\t"}{print "enhancer",$9}' > ../results/${output_dir_name}/fpkm_genes_closest_to_atac_enhancers.txt
bedtools closest -a  ../data/H1_ATAC_non_intersect_H3K27ac.bed -b ../data/H1hesc.transcript_only.tss_and_fpkm.bed | awk 'BEGIN{OFS="\t"}{print "nonenhancer",$9}' > ../results/${output_dir_name}/fpkm_genes_closest_to_atac_nonenhancers.txt
bedtools shuffle -i ../data/H1_ATAC.distal.sorted.bed -g ../data/hg19.genome -excl ../txt/excludable_ATAC_and_promoters.bed | bedtools closest -a stdin -b ../data/H1hesc.transcript_only.tss_and_fpkm.bed -d | awk 'BEGIN{OFS="\t"}{if ($7 != -1) { print "random",$9} }' > ../results/${output_dir_name}/fpkm_genes_closest_to_random_distal_regions.txt


bedtools closest -a  ../data/H1_ATAC_intersect_H3K27ac.bed -b ../data/H1hesc.transcript_only.tss_and_fpkm.bed -d | awk 'BEGIN{OFS="\t"}{if ($10 <= 10000) { print "enhancer",$9} }' > ../results/${output_dir_name}/fpkm_genes_closest_to_atac_enhancers_within10k.txt
bedtools closest -a  ../data/H1_ATAC_non_intersect_H3K27ac.bed -b ../data/H1hesc.transcript_only.tss_and_fpkm.bed -d | awk 'BEGIN{OFS="\t"}{if ($10 <= 10000) { print "nonenhancer",$9} }' > ../results/${output_dir_name}/fpkm_genes_closest_to_atac_nonenhancers_within10k.txt
bedtools shuffle -i ../data/H1_ATAC.distal.sorted.bed -g ../data/hg19.genome -excl ../txt/excludable_ATAC_and_promoters.bed | bedtools closest -a stdin -b ../data/H1hesc.transcript_only.tss_and_fpkm.bed -d | awk 'BEGIN{OFS="\t"}{if ($10 <= 10000 && $7 != -1) { print "random",$9} }' > ../results/${output_dir_name}/fpkm_genes_closest_to_random_distal_regions_within10k.txt


bedtools closest -a  ../data/H1_ATAC_intersect_H3K27ac.bed -b ../data/H1hesc.transcript_only.tss_and_fpkm.bed -d | awk 'BEGIN{OFS="\t"}{if ($10 <= 5000) { print "enhancer",$9} }' > ../results/${output_dir_name}/fpkm_genes_closest_to_atac_enhancers_within5k.txt
bedtools closest -a  ../data/H1_ATAC_non_intersect_H3K27ac.bed -b ../data/H1hesc.transcript_only.tss_and_fpkm.bed -d | awk 'BEGIN{OFS="\t"}{if ($10 <= 5000) { print "nonenhancer",$9} }' > ../results/${output_dir_name}/fpkm_genes_closest_to_atac_nonenhancers_within5k.txt
bedtools shuffle -i ../data/H1_ATAC.distal.sorted.bed -g ../data/hg19.genome -excl ../txt/excludable_ATAC_and_promoters.bed | bedtools closest -a stdin -b ../data/H1hesc.transcript_only.tss_and_fpkm.bed -d | awk 'BEGIN{OFS="\t"}{if ($10 <= 5000 && $7 != -1) { print "random",$9} }' > ../results/${output_dir_name}/fpkm_genes_closest_to_random_distal_regions_within5k.txt
Rscript expressionAnalysis.R ${output_dir_name}
## END expression analysis ###
