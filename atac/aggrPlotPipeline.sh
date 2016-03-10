bamFile=/project/umw_garberlab/tabakb/Maehr_AtacSeq/HumanES/22JUL14_PE50_C4MU1AC-B/Sample_HI_ATAC-Seq/Alignments/HI_ATAC-Seq_NoIndex_L002.Merged.cutadapt.sorted.bam
bedtools bamtobed -i $bamFile | 
awk 'BEGIN{OFS="\t"} { 
	if ($5 > 10) {
		if ($6 == "+") 
		{
			start=$2+4
			end=$3+4
		} else if ($6 == "-") {
			start=$2-5
			end=$3-5
		}
		print $1,start,end,$4,$5,$6
	}
} ' > ../data/H1_ATAC_bamToBed.filtered_and_shifted.bed

macs2 callpeak -t ../data/H1_ATAC_bamToBed.filtered_and_shifted.bed  --nomodel --shift -100 --extsize 200 -f BED -B --nomodel --SPMR -g hs -n $HOME/gnearline/ATAC_analysis/data/H1_ATAC_filtered_and_shifted
macs2 bdgcmp -t $HOME/gnearline/ATAC_analysis/data/H1_ATAC_filtered_and_shifted_treat_pileup.bdg -c $HOME/gnearline/ATAC_analysis/data/H1_ATAC_filtered_and_shifted_control_lambda.bdg -o $HOME/gnearline/ATAC_analysis/data/H1_ATAC_filtered_and_shifted.bdg -m FE
$HOME/TR/bdg2bw $HOME/gnearline/ATAC_analysis/data/H1_ATAC_filtered_and_shifted.bdg /home/si14w/TR/hg19.genome

# Do a few samples. Starting with NANOG
awk '{if ($4 == "SOX2-OCT4") print }' ../data/factorBookMotifPos.bed | bedtools intersect -a stdin -b ../data/ENCSR000AMF_H1-hESC_ChIP-seq_CTCF.gz -u > ../data/NANOG_motifsites_in_NANOG_peaks.bed
awk '{if ($4 == "SOX2-OCT4") print }' ../data/factorBookMotifPos.bed | bedtools intersect -a stdin -b ../data/ENCSR000AMF_H1-hESC_ChIP-seq_CTCF.gz -v > ../data/NANOG_motifsites_notin_NANOG_peaks.bed
bedtools slop -i ../data/NANOG_motifsites_in_NANOG_peaks.bed -b 25  -g ../data/hg19.genome | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,NR}' > ../data/NANOG_motifsites_in_NANOG_peaks_plusminus25.bed
bedtools slop -i ../data/NANOG_motifsites_notin_NANOG_peaks.bed -b 25  -g ../data/hg19.genome | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,NR}' > ../data/NANOG_motifsites_notin_NANOG_peaks_plusminus25.bed
bedtools makewindows -b ../data/NANOG_motifsites_in_NANOG_peaks_plusminus25.bed -w 1 -i srcwinnum | bigWigAverageOverBed ../data/H1_ATAC_filtered_and_shifted.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../data/ATAC_scores_NANOG_motifsites_in_NANOG_peaks_plusminus25.bed
bedtools makewindows -b ../data/NANOG_motifsites_notin_NANOG_peaks_plusminus25.bed -w 1 -i srcwinnum | bigWigAverageOverBed ../data/H1_ATAC_filtered_and_shifted.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../data/ATAC_scores_NANOG_motifsites_notin_NANOG_peaks_plusminus25.bed
awk '{ split($4,arr,"_");print arr[2]"\t"$5 }' ../data/ATAC_scores_NANOG_motifsites_in_NANOG_peaks_plusminus25.bed | sort -k1,1n | bedtools groupby -g 1 -c 2 -o mean > ../data/ATAC_profile_in_NANOGsites_in_NANOG_peaks.txt
awk '{ split($4,arr,"_");print arr[2]"\t"$5 }' ../data/ATAC_scores_NANOG_motifsites_notin_NANOG_peaks_plusminus25.bed | sort -k1,1n | bedtools groupby -g 1 -c 2 -o mean > ../data/ATAC_profile_in_NANOGsites_notin_NANOG_peaks.txt

Rscript ATAC_aggr.R ../data/ATAC_profile_in_NANOGsites_in_NANOG_peaks.txt ../data/ATAC_profile_in_NANOGsites_notin_NANOG_peaks.txt



awk '{if ($4 == "CTCF") print }' ../data/factorBookMotifPos.bed | bedtools intersect -a stdin -b ../data/ENCSR000AMF_H1-hESC_ChIP-seq_CTCF.gz -u > ../data/CTCF_motifsites_in_CTCF_peaks.bed
awk '{if ($4 == "CTCF") print }' ../data/factorBookMotifPos.bed | bedtools intersect -a stdin -b ../data/ENCSR000AMF_H1-hESC_ChIP-seq_CTCF.gz -v > ../data/CTCF_motifsites_notin_CTCF_peaks.bed
bedtools slop -i ../data/CTCF_motifsites_in_CTCF_peaks.bed -b 25  -g ../data/hg19.genome | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,NR}' > ../data/CTCF_motifsites_in_CTCF_peaks_plusminus25.bed
bedtools slop -i ../data/CTCF_motifsites_notin_CTCF_peaks.bed -b 25  -g ../data/hg19.genome | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,NR}' > ../data/CTCF_motifsites_notin_CTCF_peaks_plusminus25.bed
bedtools makewindows -b ../data/CTCF_motifsites_in_CTCF_peaks_plusminus25.bed -w 1 -i srcwinnum | bigWigAverageOverBed ../data/H1_ATAC_filtered_and_shifted.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../data/ATAC_scores_CTCF_motifsites_in_CTCF_peaks_plusminus25.bed
bedtools makewindows -b ../data/CTCF_motifsites_notin_CTCF_peaks_plusminus25.bed -w 1 -i srcwinnum | bigWigAverageOverBed ../data/H1_ATAC_filtered_and_shifted.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../data/ATAC_scores_CTCF_motifsites_notin_CTCF_peaks_plusminus25.bed
awk '{ split($4,arr,"_");print arr[2]"\t"$5 }' ../data/ATAC_scores_CTCF_motifsites_in_CTCF_peaks_plusminus25.bed | sort -k1,1n | bedtools groupby -g 1 -c 2 -o mean > ../data/ATAC_profile_in_CTCFsites_in_CTCF_peaks.txt
awk '{ split($4,arr,"_");print arr[2]"\t"$5 }' ../data/ATAC_scores_CTCF_motifsites_notin_CTCF_peaks_plusminus25.bed | sort -k1,1n | bedtools groupby -g 1 -c 2 -o mean > ../data/ATAC_profile_in_CTCFsites_notin_CTCF_peaks.txt

Rscript ATAC_aggr.R ../data/ATAC_profile_in_CTCFsites_in_CTCF_peaks.txt ../data/ATAC_profile_in_CTCFsites_notin_CTCF_peaks.txt
