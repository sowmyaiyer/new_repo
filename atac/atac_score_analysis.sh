module load IGVTools/2.3.31
igvtools tdftobedgraph /project/umw_garberlab/tabakb/Maehr_AtacSeq/HumanES/22JUL14_PE50_C4MU1AC-B/Sample_HI_ATAC-Seq/Alignments/HI_ATAC-Seq_NoIndex_L002.Merged.cutadapt.sorted.tdf ../data/H7_hESC.ATAC.bedGraph
bedGraphToBigWig ../data/H7_hESC.ATAC.bedGraph ~/TR/hg19.genome ../data/H7_hESC.ATAC.bw

cat ../data/H1_ATAC.distal.bed ../data/GencodeV19.promoter.bed > ../data/H1_ATAC_exclude_this_for_bg.bed
awk '{ print $1"\t"$2"\t"$3 }' ../data/H1_ATAC_exclude_this_for_bg.bed > tmp
mv tmp ../data/H1_ATAC_exclude_this_for_bg.bed
bedtools shuffle -i ../data/H1_ATAC.distal.bed -excl ../data/H1_ATAC_exclude_this_for_bg.bed -g ../data/hg19.genome > ../data/matched_bg_H1_ATAC.distal.bed
#awk '{ print $1"\t"$2"\t"$3"\t"$4}' ../data/H1_ATAC.distal.bed | bigWigAverageOverBed ../data/H7hESC.DNase.Stam.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/H7Hesc_DNase_Stam_signal_in_ATAC_peaks.bed
#awk '{ print $1"\t"$2"\t"$3"\t"$4}' ../data/matched_bg_H1_ATAC.distal.bed | bigWigAverageOverBed ../data/H7hESC.DNase.Stam.bw stdin ~/gnearline/RECYCLE_BIN/out2.tab -bedOut=../txt/H7Hesc_DNase_Stam_signal_in_random_regions.bed
awk '{ print $1"\t"$2"\t"$3"\t"$4}' ../data/H1_ATAC.distal.bed | bigWigAverageOverBed ../data/H7_hESC.ATAC.bw stdin ~/gnearline/RECYCLE_BIN/out.tab -bedOut=../txt/H7Hesc_ATAC_signal_in_ATAC_peaks.bed
awk '{ print $1"\t"$2"\t"$3"\t"$4}' ../data/matched_bg_H1_ATAC.distal.bed | bigWigAverageOverBed ../data/H7_hESC.ATAC.bw stdin ~/gnearline/RECYCLE_BIN/out2.tab -bedOut=../txt/H7Hesc_ATAC_signal_in_random_regions.bed

