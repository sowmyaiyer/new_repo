java -jar /share/pkg/picard/1.96/CollectInsertSizeMetrics.jar INPUT=/project/umw_garberlab/tabakb/Maehr_AtacSeq/HumanES/22JUL14_PE50_C4MU1AC-B/Sample_HI_ATAC-Seq/Alignments/HI_ATAC-Seq_NoIndex_L002.Merged.cutadapt.sorted.bam HISTOGRAM_FILE=/home/si14w/gnearline/ATAC_analysis/results/ATAC_insert_sizes.pdf OUTPUT=/home/si14w/gnearline/ATAC_analysis/results/picard_CollectInsertSizeMetrics_outFile.txt

bedtools bamtobed -i ../data/H1_ATAC_sorted_by_readname.bam.bam -bedpe |  awk '{ if ($1 != "chrM" && ($6-$2) > 38 && $8 > 10 && $1 == $4) print $0"\t"$6-$2}' > ../txt/H1_ATAC_reads_with_fragment_length.bed
bedtools intersect -a ../txt/H1_ATAC_reads_with_fragment_length.bed -b  /project/umw_garberlab/tabakb/Maehr_AtacSeq/HumanES/22JUL14_PE50_C4MU1AC-B/Sample_HI_ATAC-Seq/Alignments/macs/HI_ATAC-Seq_NoIndex_L002.Merged.cutadapt.sorted_peaks.bed -u > ../txt/ATAC_reads_in_ATAC_peaks.bed
bedtools intersect -a ../txt/ATAC_reads_in_ATAC_peaks.bed -b ../data/H1hesc.AllTFPeaks.bed -u | awk '{ print $NF"\tintersect_tf_chip"}' > ../txt/ATAC_reads_in_ATAC_peaks_intersect_tf_chip.bed
bedtools intersect -a ../txt/ATAC_reads_in_ATAC_peaks.bed -b ../data/H1hesc.AllTFPeaks.bed -v | awk '{ print $NF"\tnon_intersect_tf_chip"}' > ../txt/ATAC_reads_in_ATAC_peaks_non_intersect_tf_chip.bed

bedtools intersect -a ../txt/ATAC_reads_in_ATAC_peaks.bed -b ../data/H1_ATAC_intersect_H3K27ac.bed | awk '{ print $NF"\tenhancer"}' > ../txt/ATAC_reads_in_ATAC_peaks_intersect_enhancer.bed
bedtools intersect -a ../txt/ATAC_reads_in_ATAC_peaks.bed -b ../data/H1_ATAC_intersect_H3K27ac.bed -v | awk '{ print $NF"\tenhancer"}' > ../txt/ATAC_reads_in_ATAC_peaks_intersect_non_enhancer.bed

echo tf fragment_size > ../txt/reads_intersect_tf.txt
for tfgz in ../data/ENCS*gz
do
        tfname=`basename $tfgz | cut -d"_" -f4 | cut -d"." -f1`
        echo $tfname
	bedtools intersect -a ../txt/ATAC_reads_in_ATAC_peaks.bed -b $tfgz -u | awk -vtfname=$tfname '{ print tfname"\t"$NF }' >> ../txt/reads_intersect_tf.txt
done
