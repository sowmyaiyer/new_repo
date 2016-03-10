#sort -k5,5nr /project/umw_garberlab/tabakb/Maehr_AtacSeq/HumanES/22JUL14_PE50_C4MU1AC-B/Sample_HI_ATAC-Seq/Alignments/macs/HI_ATAC-Seq_NoIndex_L002.Merged.cutadapt.sorted_peaks.bed > ~/scratch/ATAC_sorted.bed
#rm ~/scratch/ATAC_numbers.txt
#zcat wgEncodeRegTfbsClusteredWithCellsV3.bed.gz | grep H1-hESC > H1hesc_ChIP-seq_allTFs.bed
bedtools intersect -a /project/umw_garberlab/tabakb/Maehr_AtacSeq/HumanES/22JUL14_PE50_C4MU1AC-B/Sample_HI_ATAC-Seq/Alignments/macs/HI_ATAC-Seq_NoIndex_L002.Merged.cutadapt.sorted_peaks.bed -b /project/umw_garberlab/tabakb/GSE52657_loh_lim/PeakValleyAnalyses/PeakRescoreByWindows/H3K27ac.win300.trim/trim_w50f20/H7-ESC-H3K27ac.fq.gz.vs.H7-ESC-input.fq.gz_peaks.win300.filtered.trimmed.win50fract20.bed -u | awk '{ print $NF}' > ~/TR/ATAC_scores_intersecting.txt
bedtools intersect -a /project/umw_garberlab/tabakb/Maehr_AtacSeq/HumanES/22JUL14_PE50_C4MU1AC-B/Sample_HI_ATAC-Seq/Alignments/macs/HI_ATAC-Seq_NoIndex_L002.Merged.cutadapt.sorted_peaks.bed -b /project/umw_garberlab/tabakb/GSE52657_loh_lim/PeakValleyAnalyses/PeakRescoreByWindows/H3K27ac.win300.trim/trim_w50f20/H7-ESC-H3K27ac.fq.gz.vs.H7-ESC-input.fq.gz_peaks.win300.filtered.trimmed.win50fract20.bed -v | awk '{ print $NF}' > ~/TR/ATAC_scores_non_intersecting.txt

bedtools intersect -a /project/umw_garberlab/tabakb/Maehr_AtacSeq/HumanES/22JUL14_PE50_C4MU1AC-B/Sample_HI_ATAC-Seq/Alignments/macs/HI_ATAC-Seq_NoIndex_L002.Merged.cutadapt.sorted_peaks.bed -b /project/umw_garberlab/tabakb/GSE52657_loh_lim/PeakValleyAnalyses/PeakRescoreByWindows/H3K27ac.win300.trim/trim_w50f20/H7-ESC-H3K27ac.fq.gz.vs.H7-ESC-input.fq.gz_peaks.win300.filtered.trimmed.win50fract20.bed -wao | bedtools groupby -c 4 -o distinct -full | awk '{ int_status="int"; if ($(NF-1) == 0) int_status="non_int" ; print $5"\t"int_status }' > ATAC_scores_by_intersection_status.txt

bedtools intersect -a /project/umw_garberlab/tabakb/Maehr_AtacSeq/HumanES/22JUL14_PE50_C4MU1AC-B/Sample_HI_ATAC-Seq/Alignments/macs/HI_ATAC-Seq_NoIndex_L002.Merged.cutadapt.sorted_peaks.bed -b H1hesc_ChIP-seq_allTFs.bed -u | awk '{ print $NF}' > ~/TR/ATAC_scores_intersecting_tfchip.txt
bedtools intersect -a /project/umw_garberlab/tabakb/Maehr_AtacSeq/HumanES/22JUL14_PE50_C4MU1AC-B/Sample_HI_ATAC-Seq/Alignments/macs/HI_ATAC-Seq_NoIndex_L002.Merged.cutadapt.sorted_peaks.bed -b H1hesc_ChIP-seq_allTFs.bed -v | awk '{ print $NF}' > ~/TR/ATAC_scores_non_intersecting_tfchip.txt



echo "rank_cutoff numbers" > ~/scratch/ATAC_numbers.txt
echo "rank_cutoff numbers" > ~/scratch/ATAC_tfchip_numbers.txt
echo "bin bin_desc numbers" > ~/scratch/ATAC_numbers_by_bin.txt
for cutoff in {500..150000..500}
do
	echo $cutoff
	perc=`head -$cutoff ~/scratch/ATAC_sorted.bed | bedtools intersect -a stdin -b /project/umw_garberlab/tabakb/GSE52657_loh_lim/PeakValleyAnalyses/PeakRescoreByWindows/H3K27ac.win300.trim/trim_w50f20/H7-ESC-H3K27ac.fq.gz.vs.H7-ESC-input.fq.gz_peaks.win300.filtered.trimmed.win50fract20.bed -u | wc -l`
	echo $cutoff $perc >> ~/scratch/ATAC_numbers.txt

	perc_chipseq=`head -$cutoff ~/scratch/ATAC_sorted.bed | bedtools intersect -a stdin -b H1hesc_ChIP-seq_allTFs.bed -u | wc -l`
	echo $cutoff $perc_chipseq >> ~/scratch/ATAC_tfchip_numbers.txt
done

echo "bin bin_desc numbers" > ~/scratch/ATAC_numbers_by_bin_bin1000.txt
for cutoff in {1000..150000..1000}
do
	echo $cutoff
	perc_by_bin=`head -$cutoff  ~/scratch/ATAC_sorted.bed | tail -1000 | bedtools intersect -a stdin -b /project/umw_garberlab/tabakb/GSE52657_loh_lim/PeakValleyAnalyses/PeakRescoreByWindows/H3K27ac.win300.trim/trim_w50f20/H7-ESC-H3K27ac.fq.gz.vs.H7-ESC-input.fq.gz_peaks.win300.filtered.trimmed.win50fract20.bed -u | wc -l`
        echo $((cutoff-1000))-$cutoff $perc_by_bin  >> ~/scratch/ATAC_numbers_by_bin_bin1000.txt	
done
awk '{ print (NR-1)"\t"$1"\t"$2"\t"$3 }' ~/scratch/ATAC_numbers_by_bin_bin1000.txt > ~/scratch/ATAC_numbers_by_bin_binnumbers_bin1000.txt


bedtools intersect -a /project/umw_garberlab/tabakb/Maehr_AtacSeq/HumanES/22JUL14_PE50_C4MU1AC-B/Sample_HI_ATAC-Seq/Alignments/macs/HI_ATAC-Seq_NoIndex_L002.Merged.cutadapt.sorted_peaks.bed -b H1hesc_ChIP-seq_allTFs.bed -wb > ATAC_intersect_tfChip.bed

awk '{ print $4 }' H1hesc_ChIP-seq_allTFs.bed | sort | uniq > H1hesc_tflist.txt
echo tf num_in_K27ac num_in_atac total > tf_numbers_inK27ac_and_ATAC_peaks.txt
for tf in `cut -f1 H1hesc_tflist.txt`
do
	total=`grep -c $tf H1hesc_ChIP-seq_allTFs.bed`
	num_in_K27ac=`grep $tf H1hesc_ChIP-seq_allTFs.bed | bedtools intersect -a stdin -b /project/umw_garberlab/tabakb/GSE52657_loh_lim/PeakValleyAnalyses/PeakRescoreByWindows/H3K27ac.win300.trim/trim_w50f20/H7-ESC-H3K27ac.fq.gz.vs.H7-ESC-input.fq.gz_peaks.win300.filtered.trimmed.win50fract20.bed -u | wc -l`
	num_in_atac=`grep $tf H1hesc_ChIP-seq_allTFs.bed | bedtools intersect -a stdin -b /project/umw_garberlab/tabakb/Maehr_AtacSeq/HumanES/22JUL14_PE50_C4MU1AC-B/Sample_HI_ATAC-Seq/Alignments/macs/HI_ATAC-Seq_NoIndex_L002.Merged.cutadapt.sorted_peaks.bed -u | wc -l`
	echo $tf $num_in_K27ac $num_in_atac $total >> tf_numbers_inK27ac_and_ATAC_peaks.txt
done

echo "tf number" > tfs_in_atac_peaks.txt
bedtools intersect -b H1hesc_ChIP-seq_allTFs.bed -a /project/umw_garberlab/tabakb/Maehr_AtacSeq/HumanES/22JUL14_PE50_C4MU1AC-B/Sample_HI_ATAC-Seq/Alignments/macs/HI_ATAC-Seq_NoIndex_L002.Merged.cutadapt.sorted_peaks.bed -wb | awk '{ print $9}' | sort | uniq -c | awk '{ print $2" "$1}' | sort -k2,2n >> tfs_in_atac_peaks.txt



Rscript ATAC_seq_H3K27ac_plots.R
