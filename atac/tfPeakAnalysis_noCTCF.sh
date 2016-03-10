if [[ -f ../txt/H1hesc_TF_intersection_stats_noCTCF.txt ]]; then
	rm ../txt/H1hesc_TF_intersection_stats_noCTCF.txt
fi
bedtools intersect -a ../data/H1_ATAC.distal_noCTCF.sorted.bed -b ../data/H1_and_H7_H3K27ac_merged.bed -u > ../data/H1_ATAC_noCTCF_intersect_H3K27ac.bed
bedtools intersect -a ../data/H1_ATAC.distal_noCTCF.sorted.bed -b ../data/H1_and_H7_H3K27ac_merged.bed -v > ../data/H1_ATAC_noCTCF_non_intersect_H3K27ac.bed
echo tfname peaks_intersecting_atac_enhancers peaks_intersecting_atac_nonenhancers peaks_intersecting_atac_both_enh_and_nonenh peaks_not_intersecting_atac total_peaks > ../txt/H1hesc_TF_intersection_stats_noCTCF.txt
for tfgz in ../data/ENCS*gz
do
	tfname=`basename $tfgz | cut -d"_" -f4 | cut -d"." -f1`
	echo $tfname
#	zcat $tfgz | awk -vtfname=$tfname 'BEGIN{OFS="\t"}{ print $0,tfname}' >> ../data/H1hesc.AllTFPeaks.bed
	total_peaks=`zcat $tfgz | wc -l`
	zcat $tfgz | bedtools intersect -a stdin -b ../data/H1_ATAC_noCTCF_intersect_H3K27ac.bed -wa -u > ../data/$tfname.peaks_intersecting_atac_enhancers_noCTCF.bed
	zcat $tfgz | bedtools intersect -a stdin -b ../data/H1_ATAC_noCTCF_non_intersect_H3K27ac.bed -wa -u > ../data/$tfname.peaks_intersecting_atac_nonenhancers_noCTCF.bed
	peaks_intersecting_atac_both_enh_and_nonenh=`bedtools intersect -a ../data/$tfname.peaks_intersecting_atac_enhancers_noCTCF.bed -b ../data/$tfname.peaks_intersecting_atac_nonenhancers_noCTCF.bed -f 1.0 -r | wc -l`
	peaks_intersecting_atac_enhancers=`wc -l ../data/$tfname.peaks_intersecting_atac_enhancers_noCTCF.bed | cut -d" " -f1`
	peaks_intersecting_atac_nonenhancers=`wc -l ../data/$tfname.peaks_intersecting_atac_nonenhancers_noCTCF.bed | cut -d" " -f1`
#	peaks_intersecting_atac_all=`zcat $tfgz | bedtools intersect -a stdin -b /home/si14w/gnearline/ATAC_analysis/data/H1_ATAC.distal.sorted.bed -u | wc -l`
	peaks_not_intersecting_atac=`zcat $tfgz | bedtools intersect -a stdin -b ../data/H1_ATAC.distal_noCTCF.sorted.bed -v | wc -l`
	echo $tfname $peaks_intersecting_atac_enhancers $peaks_intersecting_atac_nonenhancers $peaks_intersecting_atac_both_enh_and_nonenh $peaks_not_intersecting_atac $total_peaks >> ../txt/H1hesc_TF_intersection_stats_noCTCF.txt
done
sed 's/ /\t/g' ../txt/H1hesc_TF_intersection_stats_noCTCF.txt > ../txt/tmp
mv ../txt/tmp ../txt/H1hesc_TF_intersection_stats_noCTCF.txt


# Now asking the reverse question. How many non-enhancer ATAC peaks are bound by each TF? How many enhancer ATAC peaks are bound by each TF?

total_enhancers=`wc -l ../data/H1_ATAC.distal_noCTCF.sorted.bed | cut -d" " -f1`
total_nonenhancers=`wc -l ../data/H1_ATAC_noCTCF_non_intersect_H3K27ac.bed | cut -d" " -f1`
echo  TF enhancers_bound_to_tf total_enhancers totalTFPeaks > ../txt/enhancer_intersecting_tf_noCTCF.txt
echo  TF nonenhancers_bound_to_tf total_nonenhancers totalTFPeaks > ../txt/nonenhancer_intersecting_tf_noCTCF.txt
for tfgz in ../data/ENCS*gz
do
	tfname=`basename $tfgz | cut -d"_" -f4 | cut -d"." -f1`
	totalPeaks=`zcat $tfgz| wc -l | cut -d" " -f1`
	enhancers_intersecting_tf=`zcat $tfgz | bedtools intersect -a ../data/H1_ATAC_noCTCF_intersect_H3K27ac.bed -b stdin -u | wc -l`
	nonenhancers_intersecting_tf=`zcat $tfgz | bedtools intersect -a ../data/H1_ATAC_noCTCF_non_intersect_H3K27ac.bed -b stdin -u | wc -l`
	echo $tfname $enhancers_intersecting_tf $total_enhancers $totalPeaks >> ../txt/enhancer_intersecting_tf_noCTCF.txt
	echo $tfname $nonenhancers_intersecting_tf $total_nonenhancers $totalPeaks >> ../txt/nonenhancer_intersecting_tf_noCTCF.txt
done

bedtools intersect -a ../data/H1_ATAC.distal_noCTCF.sorted.bed -b ../data/H1hesc.AllTFPeaks.bed -wa -wb | awk 'BEGIN{OFS="\t"}{ print $16,$5}' > ../txt/TFs_and_ATAC_scores_noCTCF.txt

bedtools intersect -a ../data/H1hesc.AllTFPeaks.bed -b ../data/H1_ATAC_noCTCF_intersect_H3K27ac.bed | awk 'BEGIN{OFS="\t"}{ print $11,$7,"enhancer_ATAC"}' > ../txt/H1_TFPeaks_in_enhancer_atac_peaks_noCTCF.txt
bedtools intersect -a ../data/H1hesc.AllTFPeaks.bed -b ../data/H1_ATAC_noCTCF_non_intersect_H3K27ac.bed | awk 'BEGIN{OFS="\t"}{ print $11,$7,"non_enhancer_ATAC" }' > ../txt/H1_TFPeaks_in_nonenhancer_atac_peaks_noCTCF.txt
bedtools intersect -a ../data/H1hesc.AllTFPeaks.bed -b ../data/H1_ATAC.distal_noCTCF.sorted.bed -v | awk 'BEGIN{OFS="\t"}{ print $11,$7,"non_ATAC"}' > ../txt/H1_TFPeaks_not_in_atac_peaks_noCTCF.txt
