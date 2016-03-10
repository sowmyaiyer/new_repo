for time in {"6h","12h","24h","48h","mock12h"}
do
	echo """
	module load bowtie2/2-2.1.0
	module load tophat/2.0.14
	mkdir  /home/si14w/gnearline/flu/tophat_out_${time}
	tophat -p 8 --library-type fr-firststrand -o /home/si14w/gnearline/flu/tophat_out_${time} /home/si14w/gnearline/flu/bowtie_out/hg19_and_fluWithSnp_combined /farline/umw_flustore_schiffer/rnaseq/04SEP15.A549_B59_${time}.fastq.gz 
	""" > ../bsubFiles/runtophat_rnaseq.${time}.bsub
done
