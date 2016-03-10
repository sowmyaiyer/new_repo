while read line
do
	r1_file=`echo $line | awk '{ print $2}'`
	r2_file=`echo $line | awk '{ print $3}'`
	sample=`echo $line | awk '{ print $1}'`
	cat ~/gnearline/hdc/chip1124/${r1_file} >> ~/gnearline/hdc/chip1124/${sample}_R1.gz
	cat ~/gnearline/hdc/chip1124/${r2_file} >> ~/gnearline/hdc/chip1124/${sample}_R2.gz
	ls -ltr ~/gnearline/hdc/chip1124/${sample}_R[1,2].gz
done<../txt/chip1124_sampleInfo.txt

module load fastqc/0.10.1
for  fastq_gz in  ~/gnearline/hdc/chip1124/*R[1,2].gz
do
	 fastqc $fastq_gz -o /home/si14w/gnearline/hdc/chipseq_fastqc
done
