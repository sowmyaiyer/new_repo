for time in {"A","B","C","D","E","F","G","H"}
do
	echo """Rscript findTSS.R home/si14w/gnearline/flu/bowtie_out/112213flucap_${time}.sorted.bam ctss_${time}""" > ../bsubFiles/findTss_${time}.bsub
done
