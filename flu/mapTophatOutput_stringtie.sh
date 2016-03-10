for time in {"6h","12h","24h","48h","mock12h"}
do
	echo """
	module load stringtie/1.1.2
	mkdir /home/si14w/gnearline/flu/ballgown_input_${time}
	stringtie -p 8 -b /home/si14w/gnearline/flu/ballgown_input_${time} -o /home/si14w/gnearline/flu/stringtie_out/stringtie_out_${time} /home/si14w/gnearline/flu/tophat_out_${time}/accepted_hits.bam 
	""" > ../bsubFiles/runstringtie_rnaseq.${time}.bsub
done
