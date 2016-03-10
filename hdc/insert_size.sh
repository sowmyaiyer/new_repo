for ab in {"H3K4me3","H3K27ac","H3K4me1","IgG","input"}
do
	time="0h"
	java -jar /share/pkg/picard/1.96/CollectInsertSizeMetrics.jar INPUT=/home/si14w/gnearline/hdc/bwa_out/E37_6d_${time}_${ab}_1mln_cells.pe.sorted.bam HISTOGRAM_FILE=../results/chip0929_qc_insert_size_${time}_${ab}.pdf OUTPUT=../results/chip0929_qc_insert_size_${time}_${ab}.txt
	echo done $ab
done
