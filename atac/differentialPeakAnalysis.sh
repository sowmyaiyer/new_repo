#10/5
sort -k1,1 ~/gnearline/hdc/txt/D50_6d_LPS.tsv > ../txt/D50_6d_LPS.sorted_by_genename.tsv
bedtools intersect -a ../hdc_peaks/merged_6d_4h_LPS_peaks.bed -b ../hdc_peaks/merged_6d_1h_LPS_peaks.bed -v | bedtools intersect -a stdin -b ../hdc_peaks/merged_6d_0h_LPS_peaks.bed -v | bedtools intersect -a stdin -b ~/gnearline/common_scripts/hg19_refseq_rsem_corrected_tss.sorted.5K.bed -wb | sort -k9,9 > ../txt/new_peaks_4h_LPS.bed
join ../txt/new_peaks_4h_LPS.bed ../txt/D50_6d_LPS.sorted_by_genename.tsv -1 9 -2 1 | sed 's/ /\t/g' | sort -k6,6nr -k14,14nr


bedtools multiinter -i ../hdc_peaks/merged_6d_0h_LPS_peaks.bed ../hdc_peaks/merged_6d_1h_LPS_peaks.bed ../hdc_peaks/merged_6d_4h_LPS_peaks.bed -names 0h 1h 4h | grep "0h,1h,4h" | bedtools intersect -a  stdin -b ~/gnearline/common_scripts/hg19_refseq_rsem_corrected_tss.sorted.5K.bed -wb | sort -k12,12 > ../txt/peaks_all_3_time_points.bed
join ../txt/peaks_all_3_time_points.bed ../txt/D50_6d_LPS.sorted_by_genename.tsv -1 12 -2 1 | sed 's/ /\t/g' | awk '{ print $0"\t"($16+1)/($15+1) }' | sort -k21,21nr | more
bedtools multiinter -i ../hdc_peaks/merged_6d_0h_LPS_peaks.bed ../hdc_peaks/merged_6d_1h_LPS_peaks.bed ../hdc_peaks/merged_6d_4h_LPS_peaks.bed -names 0h 1h 4h | awk '{ if ($5 == "0h,1h") print }'  | bedtools intersect -a  stdin -b ~/gnearline/common_scripts/hg19_refseq_rsem_corrected_tss.sorted.5K.bed -wb | sort -k12,12 > ../txt/peaks_0h_1h.bed
join ../txt/peaks_0h_1h.bed ../txt/D50_6d_LPS.sorted_by_genename.tsv -1 12 -2 1 | sed 's/ /\t/g' | awk '{ print $0"\t"($16+1)/($15+1) }' | sort -k21,21nr | more

#10/6
bedtools window -a ../hdc_peaks/merged_6d_1h_LPS_peaks.bed -w 1000 -b ../hdc_peaks/merged_6d_0h_LPS_peaks.bed -v | bedtools intersect -a stdin -b  ~/gnearline/common_scripts/hg19_refseq_rsem_corrected_tss.sorted.5K.bed -wb | sort -k9,9 > ../txt/peaks_1h_new.bed
join ../txt/peaks_1h_new.bed   ../txt/D50_6d_LPS.sorted_by_genename.tsv -1 9 -2 1 | sed 's/ /\t/g' | awk '{ print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$1"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17}' | sort -k5,5nr > ../txt/peaks_1h_new_intersect_genes.txt
