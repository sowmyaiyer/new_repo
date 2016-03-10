while read sample
do
#	cat /farline/umw_flustore_schiffer/eIF4E_pulldown_RNAseq/fastq_untrimmed/${sample}_L00*R1_001.fastq.gz	> /farline/umw_flustore_schiffer/eIF4E_pulldown_RNAseq/fastq_untrimmed/${sample}.R1.fastq.gz
echo """
	module load fastqc/0.10.1
        fastqc  /farline/umw_flustore_schiffer/eIF4E_pulldown_RNAseq/fastq_untrimmed/${sample}.R1.fastq.gz -o /home/si14w/gnearline/flu/fastqc_results/
""" > ../bsubFiles/runFastqc.${sample}.bsub
done<samples.txt
