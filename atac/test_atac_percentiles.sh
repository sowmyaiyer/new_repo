awk '{ print $NF}' /project/umw_garberlab/tabakb/Maehr_AtacSeq/HumanES/22JUL14_PE50_C4MU1AC-B/Sample_HI_ATAC-Seq/Alignments/macs/HI_ATAC-Seq_NoIndex_L002.Merged.cutadapt.sorted_peaks.bed > ATAC_scores_all.txt
Rscript getAtacScorePercentiles.R /home/si14w/scratch/ATAC_scores_all.txt /home/si14w/scratch/atac_quantile_values.txt
echo "bin numbers total" > atac_scores_percentiles.txt
for n in {1..100}
do
	echo $n
	greater_than=`head -$((n+1)) ~/scratch/atac_quantile_values.txt | tail -1`
	less_than=`head -$n ~/scratch/atac_quantile_values.txt | tail -1`
	echo $greater_than $less_than 
	awk -vgreater_than=$greater_than -vless_than=$less_than '{ if ($NF > greater_than && $NF < less_than) print }' /project/umw_garberlab/tabakb/Maehr_AtacSeq/HumanES/22JUL14_PE50_C4MU1AC-B/Sample_HI_ATAC-Seq/Alignments/macs/HI_ATAC-Seq_NoIndex_L002.Merged.cutadapt.sorted_peaks.bed | wc -l
	total=`awk \-vgreater_than=$greater_than \-vless_than=$less_than '{ if ($NF > greater_than && $NF < less_than) print }' /project/umw_garberlab/tabakb/Maehr_AtacSeq/HumanES/22JUL14_PE50_C4MU1AC-B/Sample_HI_ATAC-Seq/Alignments/macs/HI_ATAC-Seq_NoIndex_L002.Merged.cutadapt.sorted_peaks.bed | wc -l`
	number_intersecting=`awk \-vgreater_than=$greater_than \-vless_than=$less_than '{ if ($NF > greater_than && $NF < less_than) print }' /project/umw_garberlab/tabakb/Maehr_AtacSeq/HumanES/22JUL14_PE50_C4MU1AC-B/Sample_HI_ATAC-Seq/Alignments/macs/HI_ATAC-Seq_NoIndex_L002.Merged.cutadapt.sorted_peaks.bed | bedtools intersect -a stdin -b /project/umw_garberlab/tabakb/GSE52657_loh_lim/PeakValleyAnalyses/PeakRescoreByWindows/H3K27ac.win300.trim/trim_w50f20/H7-ESC-H3K27ac.fq.gz.vs.H7-ESC-input.fq.gz_peaks.win300.filtered.trimmed.win50fract20.bed -u | wc -l`
	echo $greater_than $number_intersecting $total >> atac_scores_percentiles.txt
done

bedtools slop -i /home/si14w/nearline/enhancer_predictions/gencodeV19_TSS.bed -l 1000 -r 500 -s -g ~/TR/hg19.genome > GencodeV19.promoter.bed

bedtools intersect -a /project/umw_garberlab/tabakb/Maehr_AtacSeq/HumanES/22JUL14_PE50_C4MU1AC-B/Sample_HI_ATAC-Seq/Alignments/macs/HI_ATAC-Seq_NoIndex_L002.Merged.cutadapt.sorted_peaks.bed -b GencodeV19.promoter.bed -u | bedtools intersect -a stdin -b /project/umw_garberlab/tabakb/GSE52657_loh_lim/PeakValleyAnalyses/PeakRescoreByWindows/H3K27ac.win300.trim/trim_w50f20/H7-ESC-H3K27ac.fq.gz.vs.H7-ESC-input.fq.gz_peaks.win300.filtered.trimmed.win50fract20.bed -u | awk '{ print $NF}' > ~/TR/ATAC_scores_intersecting_H3K27ac.proximal.txt
bedtools intersect -a /project/umw_garberlab/tabakb/Maehr_AtacSeq/HumanES/22JUL14_PE50_C4MU1AC-B/Sample_HI_ATAC-Seq/Alignments/macs/HI_ATAC-Seq_NoIndex_L002.Merged.cutadapt.sorted_peaks.bed -b GencodeV19.promoter.bed -v | bedtools intersect -a stdin -b /project/umw_garberlab/tabakb/GSE52657_loh_lim/PeakValleyAnalyses/PeakRescoreByWindows/H3K27ac.win300.trim/trim_w50f20/H7-ESC-H3K27ac.fq.gz.vs.H7-ESC-input.fq.gz_peaks.win300.filtered.trimmed.win50fract20.bed -u | awk '{ print $NF}' > ~/TR/ATAC_scores_intersecting_H3K27ac.distal.txt

bedtools intersect -a /project/umw_garberlab/tabakb/Maehr_AtacSeq/HumanES/22JUL14_PE50_C4MU1AC-B/Sample_HI_ATAC-Seq/Alignments/macs/HI_ATAC-Seq_NoIndex_L002.Merged.cutadapt.sorted_peaks.bed -b GencodeV19.promoter.bed -u | bedtools intersect -a stdin -b /project/umw_garberlab/tabakb/GSE52657_loh_lim/PeakValleyAnalyses/PeakRescoreByWindows/H3K27ac.win300.trim/trim_w50f20/H7-ESC-H3K27ac.fq.gz.vs.H7-ESC-input.fq.gz_peaks.win300.filtered.trimmed.win50fract20.bed -v | awk '{ print $NF}' > ~/TR/ATAC_scores_nonintersecting_H3K27ac.proximal.txt
bedtools intersect -a /project/umw_garberlab/tabakb/Maehr_AtacSeq/HumanES/22JUL14_PE50_C4MU1AC-B/Sample_HI_ATAC-Seq/Alignments/macs/HI_ATAC-Seq_NoIndex_L002.Merged.cutadapt.sorted_peaks.bed -b GencodeV19.promoter.bed -v | bedtools intersect -a stdin -b /project/umw_garberlab/tabakb/GSE52657_loh_lim/PeakValleyAnalyses/PeakRescoreByWindows/H3K27ac.win300.trim/trim_w50f20/H7-ESC-H3K27ac.fq.gz.vs.H7-ESC-input.fq.gz_peaks.win300.filtered.trimmed.win50fract20.bed -v | awk '{ print $NF}' > ~/TR/ATAC_scores_nonintersecting_H3K27ac.distal.txt

Rscript ATAC_scores_proximal_distal.R
