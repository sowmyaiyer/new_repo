for file in `ls ~/gnearline/hdc/chip0929/C*/*gz`
do
	mv $file ~/gnearline/hdc/chip0929/
done
while read line
do
	r1_file=`echo $line | awk '{ print $2}'`
	r2_file=`echo $line | awk '{ print $3}'`
	sample=`echo $line | awk '{ print $1}'`
	cat ~/gnearline/hdc/chip0929/${r1_file} >> ~/gnearline/hdc/chip0929/${sample}_R1.gz
	cat ~/gnearline/hdc/chip0929/${r2_file} >> ~/gnearline/hdc/chip0929/${sample}_R2.gz
	ls -ltr ~/gnearline/hdc/chip0929/${sample}_R[1,2].gz
done<../txt/chip0929_sampleInfo.txt

module load fastqc/0.10.1
for  fastq_gz in  ~/gnearline/hdc/chip0929/*R[1,2].gz
do
	 fastqc $fastq_gz -o /home/si14w/gnearline/hdc/chipseq_fastqc
done
