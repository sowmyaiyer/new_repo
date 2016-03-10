module load fastqc/0.10.1
while read line
do
	barcode=`echo $line | cut -d" " -f2`
	sample=`echo $line | cut -d" " -f1`
	cat  /farline/umw_flustore_schiffer/unknown_maybeRNAseq/24SEP15_PE100_C7DMPAC/Sample_B59-A549mix/B59-A549mix_${barcode}_*_R1_*fastq.gz > /farline/umw_flustore_schiffer/unknown_maybeRNAseq/24SEP15.${sample}.R1.fastq.gz
	cat  /farline/umw_flustore_schiffer/unknown_maybeRNAseq/24SEP15_PE100_C7DMPAC/Sample_B59-A549mix/B59-A549mix_${barcode}_*_R2_*fastq.gz > /farline/umw_flustore_schiffer/unknown_maybeRNAseq/24SEP15.${sample}.R2.fastq.gz
	fastqc /farline/umw_flustore_schiffer/unknown_maybeRNAseq/24SEP15.${sample}.R1.fastq.gz -o /home/si14w/gnearline/flu/rnaseq_fastqc_newpe
	fastqc /farline/umw_flustore_schiffer/unknown_maybeRNAseq/24SEP15.${sample}.R2.fastq.gz -o /home/si14w/gnearline/flu/rnaseq_fastqc_newpe
	echo done $sample
done</home/si14w/gnearline/flu/scripts/flu_rnaseq_barcodes.txt
