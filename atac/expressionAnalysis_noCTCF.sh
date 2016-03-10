bedtools intersect -a ../data/H1_ATAC.distal_noCTCF.sorted.bed -b ../data/H1_and_H7_H3K27ac_merged.bed -u > ../data/H1_ATAC_noCTCF_intersect_H3K27ac.bed
bedtools intersect -a ../data/H1_ATAC.distal_noCTCF.sorted.bed -b ../data/H1_and_H7_H3K27ac_merged.bed -v > ../data/H1_ATAC_noCTCF_non_intersect_H3K27ac.bed
cat ../data/GencodeV19.promoter.bed ../data/H1_ATAC.distal.sorted.bed | awk 'BEGIN{OFS="\t"}{ print $1,$2,$3}' > ../txt/random_regions_excluding_ATAC_and_promoters.bed

awk '{ if ($3 == "transcript") print }' ../data/H1hesc.transcript.quant.gtf  | awk -f extractFpkmWithTss.awk > ../data/H1hesc.transcript_only.tss_and_fpkm.bed

bedtools closest -a  ../data/H1_ATAC_noCTCF_intersect_H3K27ac.bed -b ../data/H1hesc.transcript_only.tss_and_fpkm.bed | awk 'BEGIN{OFS="\t"}{print "enhancer",$9}' > ../txt/fpkm_genes_closest_to_atac_enhancers_noCTCF.txt
bedtools closest -a  ../data/H1_ATAC_noCTCF_non_intersect_H3K27ac.bed -b ../data/H1hesc.transcript_only.tss_and_fpkm.bed | awk 'BEGIN{OFS="\t"}{print "nonenhancer",$9}' > ../txt/fpkm_genes_closest_to_atac_nonenhancers_noCTCF.txt
bedtools shuffle -i /home/si14w/gnearline/ATAC_analysis/data/H1_ATAC.distal_noCTCF.sorted.bed -g ../data/hg19.genome -excl ../txt/random_regions_excluding_ATAC_and_promoters.bed | bedtools closest -a stdin -b ../data/H1hesc.transcript_only.tss_and_fpkm.bed -d | awk 'BEGIN{OFS="\t"}{if ($7 != -1) { print "random",$9} }' > ../txt/fpkm_genes_closest_to_random_distal_regions_noCTCF.txt


bedtools slop -i ../data/H1_ATAC_noCTCF_intersect_H3K27ac.bed -b 5000 -g ../data/hg19.genome | bedtools intersect -a stdin -b ../data/H1hesc.transcript_only.tss_and_fpkm.bed -wa -wb | awk  'BEGIN{OFS="\t"}{print "enhancer",$9}' > ../txt/enhancer_fpkm_all_genes_within_5000_noCTCF.txt 
bedtools slop -i ../data/H1_ATAC_noCTCF_non_intersect_H3K27ac.bed -b 5000 -g ../data/hg19.genome | bedtools intersect -a stdin -b ../data/H1hesc.transcript_only.tss_and_fpkm.bed -wa -wb | awk 'BEGIN{OFS="\t"}{print "nonenhancer",$9}' > ../txt/nonenhancer_fpkm_all_genes_within_5000_noCTCF.txt

#bedtools slop -i ../data/H1_ATAC_noCTCF_intersect_H3K27ac.bed -b 5000 -g ../data/hg19.genome | bedtools intersect -b stdin -a ../data/H1hesc.transcript_only.tss_and_fpkm.bed -wa | awk -OFS="\t" '{print $1,$2,$3,$4 }' > ../txt/gene_tss_and_fpkm_within5k_of_enhancers.bed
#bedtools slop -i ../data/H1_ATAC_noCTCF_non_intersect_H3K27ac.bed -b 5000 -g ../data/hg19.genome | bedtools intersect -b stdin -a ../data/H1hesc.transcript_only.tss_and_fpkm.bed -wa | awk -OFS="\t" '{print $1,$2,$3,$4 }' > ../txt/gene_tss_and_fpkm_within5k_of_nonenhancers.bed

bedtools closest -a  ../data/H1_ATAC_noCTCF_intersect_H3K27ac.bed -b ../data/H1hesc.transcript_only.tss_and_fpkm.bed -d | awk 'BEGIN{OFS="\t"}{if ($10 <= 10000) { print "enhancer",$9} }' > ../txt/fpkm_genes_closest_to_atac_enhancers_within10k_noCTCF.txt
bedtools closest -a  ../data/H1_ATAC_noCTCF_non_intersect_H3K27ac.bed -b ../data/H1hesc.transcript_only.tss_and_fpkm.bed -d | awk 'BEGIN{OFS="\t"}{if ($10 <= 10000) { print "nonenhancer",$9} }' > ../txt/fpkm_genes_closest_to_atac_nonenhancers_within10k_noCTCF.txt
bedtools shuffle -i /home/si14w/gnearline/ATAC_analysis/data/H1_ATAC.distal_noCTCF.sorted.bed -g ../data/hg19.genome -excl ../txt/random_regions_excluding_ATAC_and_promoters.bed | bedtools closest -a stdin -b ../data/H1hesc.transcript_only.tss_and_fpkm.bed -d | awk 'BEGIN{OFS="\t"}{if ($10 <= 10000 && $7 != -1) { print "random",$9} }' > ../txt/fpkm_genes_closest_to_random_distal_regions_within10k_noCTCF.txt


bedtools closest -a  ../data/H1_ATAC_noCTCF_intersect_H3K27ac.bed -b ../data/H1hesc.transcript_only.tss_and_fpkm.bed -d | awk 'BEGIN{OFS="\t"}{if ($10 <= 5000) { print "enhancer",$9} }' > ../txt/fpkm_genes_closest_to_atac_enhancers_within5k_noCTCF.txt
bedtools closest -a  ../data/H1_ATAC_noCTCF_non_intersect_H3K27ac.bed -b ../data/H1hesc.transcript_only.tss_and_fpkm.bed -d | awk 'BEGIN{OFS="\t"}{if ($10 <= 5000) { print "nonenhancer",$9} }' > ../txt/fpkm_genes_closest_to_atac_nonenhancers_within5k_noCTCF.txt
bedtools shuffle -i /home/si14w/gnearline/ATAC_analysis/data/H1_ATAC.distal_noCTCF.sorted.bed -g ../data/hg19.genome -excl ../txt/random_regions_excluding_ATAC_and_promoters.bed | bedtools closest -a stdin -b ../data/H1hesc.transcript_only.tss_and_fpkm.bed -d | awk 'BEGIN{OFS="\t"}{if ($10 <= 5000 && $7 != -1) { print "random",$9} }' > ../txt/fpkm_genes_closest_to_random_distal_regions_within5k_noCTCF.txt


bedtools closest -a  ../data/H1_ATAC_noCTCF_intersect_H3K27ac.bed -b ../data/H1hesc.transcript_only.tss_and_fpkm.bed -d > ../txt/enhancers_and_closest_genes_noCTCF.bed
bedtools closest -a  ../data/H1_ATAC_noCTCF_non_intersect_H3K27ac.bed -b ../data/H1hesc.transcript_only.tss_and_fpkm.bed -d > ../txt/nonenhancers_and_closest_genes_noCTCF.bed

bedtools shuffle -i /home/si14w/gnearline/ATAC_analysis/data/H1_ATAC.distal_noCTCF.sorted.bed -g ../data/hg19.genome -excl ../data/GencodeV19.promoter.bed | bedtools closest -a stdin -b ../data/H1hesc.transcript_only.tss_and_fpkm.bed -d > ../txt/fpkm_genes_closest_to_random_distal_regions_noCTCF.bed
