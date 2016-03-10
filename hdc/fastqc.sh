for sample in  `awk '{ print $1}' ../txt/chipnew_sampleInfo.txt | sort | uniq`
do
	echo """
        module load fastqc/0.10.1
        fastqc  ~/gnearline/hdc/chipnew/${sample}_R1.gz -o /home/si14w/gnearline/hdc/chipseq_fastqc
        fastqc  ~/gnearline/hdc/chipnew/${sample}_R2.gz -o /home/si14w/gnearline/hdc/chipseq_fastqc
""" > ../bsubFiles/fastqc_${sample}.bsub
done
