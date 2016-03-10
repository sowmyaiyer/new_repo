module load fastqc/0.10.1
while read line
do
	barcode=`echo $line | cut -d" " -f2`
	sample=`echo $line | cut -d" " -f1`
	cat /farline/umw_flustore_schiffer/rnaseq/04SEP15_SR100_C7RMHAC/Sample_B59-A549_mix/B59-A549_mix_${barcode}*fastq.gz > /farline/umw_flustore_schiffer/rnaseq/04SEP15.${sample}.fastq.gz
	fastqc /farline/umw_flustore_schiffer/rnaseq/04SEP15.${sample}.fastq.gz -o /home/si14w/gnearline/flu/rnaseq_fastqc
	echo done $sample
done</home/si14w/gnearline/flu/scripts/flu_rnaseq_barcodes.txt
