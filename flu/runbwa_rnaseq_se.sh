#bwa index -p ../bwa_out/hg19_and_Brisbane59_combined.fa ../hg19_and_Brisbane59_combined.fa
for time in {"6h","12h","24h","48h","mock12h"}
do
	echo """
	#bwa aln -q 20 -t 8 ../bwa_out/hg19_and_Brisbane59_combined.fa /farline/umw_flustore_schiffer/rnaseq/04SEP15.A549_B59_${time}.fastq.gz > ../bwa_out/04SEP15.A549_B59_${time}.sai
	#bwa samse ../bwa_out/hg19_and_Brisbane59_combined.fa ../bwa_out/04SEP15.A549_B59_${time}.sai /farline/umw_flustore_schiffer/rnaseq/04SEP15.A549_B59_${time}.fastq.gz > ../bwa_out/04SEP15.A549_B59_${time}.sam
	samtools view -bS -o ../bwa_out/04SEP15.A549_B59_${time}.bam ../bwa_out/04SEP15.A549_B59_${time}.sam
	samtools sort -T ../bwa_out/04SEP15.A549_B59_${time}.sorted -o ../bwa_out/04SEP15.A549_B59_${time}.sorted.bam ../bwa_out/04SEP15.A549_B59_${time}.bam
#	samtools index ../bwa_out/04SEP15.A549_B59_${time}.sorted.bam
	java -jar /share/pkg/picard/1.96/MarkDuplicates.jar INPUT=04SEP15.A549_B59_${time}.sorted.bam OUTPUT=04SEP15.A549_B59_${time}.sorted.deduped.bam METRICS_FILE=04SEP15.A549_B59_${time}.sorted.bam.duplicates.metrics.out REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT
	samtools index 04SEP15.A549_B59_${time}.sorted.deduped.bam
	rm  ../bwa_out/04SEP15.A549_B59_${time}.sai ../bwa_out/04SEP15.A549_B59_${time}.sam
	""" > ../bsubFiles/flu_rnaseq_bwa_${time}.bsub
done
