for file in `ls ~/gnearline/hdc/chipnew/ChIP*/*gz`
do
	mv $file ~/gnearline/hdc/chipnew/
done
while read line
do
	r1_file=`echo $line | awk '{ print $2}'`
	r2_file=`echo $line | awk '{ print $3}'`
	sample=`echo $line | awk '{ print $1}'`
	cat ~/gnearline/hdc/chipnew/${r1_file} >> ~/gnearline/hdc/chipnew/${sample}_R1.gz
	cat ~/gnearline/hdc/chipnew/${r2_file} >> ~/gnearline/hdc/chipnew/${sample}_R2.gz
	ls -ltr ~/gnearline/hdc/chipnew/${sample}_R[1,2].gz
done<../txt/chipnew_sampleInfo.txt

module load fastqc/0.10.1
for  fastq_gz in  ~/gnearline/hdc/chipnew/*R[1,2].gz
do
	 fastqc $fastq_gz -o /home/si14w/gnearline/hdc/chipseq_fastqc
done
