#ls /farline/umw_flustore_schiffer/hdc_chip/Chip_retest/Chip_retest-28848823/*/*fastq.gz | awk -F"/" '{ print $NF}'  | cut -d"/" -f3 | cut -d"_" -f1,2 | sort | uniq > /home/si14w/gnearline/hdc/txt/hdc_patientsamples_info_02_26.txt
for sample in `awk '{ print $1}' /home/si14w/gnearline/hdc/txt/hdc_patientsamples_info_03_05.txt`
do
	echo $sample 
        cat /farline/umw_flustore_schiffer/hdc_chip/Chip_retest_2/ChIP_seq-29006978/*/${sample}_*R1_001.fastq.gz > /farline/umw_garberlab/human/DC/ChIP-Seq/patient_samples/chip_retest_2/${sample}.R1.fastq.gz
        cat /farline/umw_flustore_schiffer/hdc_chip/Chip_retest_2/ChIP_seq-29006978/*/${sample}_*R2_001.fastq.gz > /farline/umw_garberlab/human/DC/ChIP-Seq/patient_samples/chip_retest_2/${sample}.R2.fastq.gz
	ls /farline/umw_garberlab/human/DC/ChIP-Seq/patient_samples/chip_retest_2/${sample}.R1.fastq.gz /farline/umw_garberlab/human/DC/ChIP-Seq/patient_samples/chip_retest_2/${sample}.R2.fastq.gz
	echo """
	module load fastqc/0.10.1
	fastqc /farline/umw_garberlab/human/DC/ChIP-Seq/patient_samples/chip_retest_2/${sample}.R1.fastq.gz -o /home/si14w/gnearline/hdc/chipseq_fastqc/chip_retest2/
	fastqc /farline/umw_garberlab/human/DC/ChIP-Seq/patient_samples/chip_retest_2/${sample}.R2.fastq.gz -o /home/si14w/gnearline/hdc/chipseq_fastqc/chip_retest2/
	""" > ../bsubFiles/fastqc_chipretest2_$sample.bsub
done
