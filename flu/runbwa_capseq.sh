#bwa index -p ../bwa_out/hg19_and_Brisbane59_combined.fa ../hg19_and_Brisbane59_combined.fa
for time in {"A","B","C","D","F"}
do
	echo """
	bwa aln -t 8 ../bwa_out/hg19_and_Brisbane59_combined.fa /project/umw_garberlab/narayanan/112213flucap/rawData_toFluAndtoHuman/${time}_Nocodes.fastq.gz > ../bwa_out/112213flucap_${time}.sai
	bwa samse ../bwa_out/hg19_and_Brisbane59_combined.fa ../bwa_out/112213flucap_${time}.sai /project/umw_garberlab/narayanan/112213flucap/rawData_toFluAndtoHuman/${time}_Nocodes.fastq.gz > ../bwa_out/112213flucap_${time}.sam
	samtools view -bS -o ../bwa_out/112213flucap_${time}.bam ../bwa_out/112213flucap_${time}.sam
	samtools sort -T ../bwa_out/112213flucap_${time}.sorted -o ../bwa_out/112213flucap_${time}.sorted.bam ../bwa_out/112213flucap_${time}.bam
	samtools index ../bwa_out/112213flucap_${time}.sorted.bam
	#java -jar /share/pkg/picard/1.96/MarkDuplicates.jar INPUT=112213flucap_${time}.sorted.bam OUTPUT=112213flucap_${time}.sorted.deduped.bam METRICS_FILE=112213flucap_${time}.sorted.bam.duplicates.metrics.out REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT
	#samtools index 112213flucap_${time}.sorted.deduped.bam
	rm  ../bwa_out/112213flucap_${time}.sai ../bwa_out/112213flucap_${time}.sam
	""" > ../bsubFiles/flu_capseq_bwa_${time}.bsub
done
