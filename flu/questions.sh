echo segment reads_aligning_to_sense_strand_first_150 reads_aligning_to_antisense_strand_last_150 total_reads_aligning_to_sense_strand	total_reads_aligning_to_antisense_strand > flu_segmentwise_alignment_info.txt
for segment in {1..8}
do	
	segname="seg"$segment"_B59FULL"
	seglen=`awk -vseg=$segment '{ if ($1 == seg) print $2}' segment_lengths.txt`
	echo $segname 
	reads_aligning_to_antisense_strand_last_150=`samtools view -f 16 ../data/112213_barcodeC-sorted.bam | awk -vlen=$seglen -vname=$segname '{ if (($3 == name) && ($4 > len-46-150)) print }' | wc -l`
	reads_aligning_to_sense_strand_first_150=`samtools view -F 16 ../data/112213_barcodeC-sorted.bam | awk -vlen=$seglen -vname=$segname '{ if (($3 == name) && ($4 < 150)) print }' | wc -l`
	total_reads_aligning_to_sense_strand=`samtools view -F 16 ../data/112213_barcodeC-sorted.bam | awk -vname=$segname '{ if ($3 == name) print }' | wc -l`
	total_reads_aligning_to_antisense_strand=`samtools view -f 16 ../data/112213_barcodeC-sorted.bam | awk -vname=$segname '{ if ($3 == name) print }' |  wc -l`	
	echo $segname $reads_aligning_to_sense_strand_first_150 $reads_aligning_to_antisense_strand_last_150 $total_reads_aligning_to_sense_strand $total_reads_aligning_to_antisense_strand >> flu_segmentwise_alignment_info.txt
done
