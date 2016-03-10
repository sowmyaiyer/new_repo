for time in {"A","B","C","D","E","F","G","H"}
do
	echo """
	samtools view -b -q 10 -o  /home/si14w/gnearline/flu/bowtie_out/112213flucap_${time}.filtered.bam /home/si14w/gnearline/flu/bowtie_out/112213flucap_${time}.sorted.bam 
	samtools sort -T ../bowtie_out/112213flucap_${time}.filtered.sorted -o ../bowtie_out/112213flucap_${time}.filtered.sorted.bam /home/si14w/gnearline/flu/bowtie_out/112213flucap_${time}.filtered.bam
        samtools index ../bowtie_out/112213flucap_${time}.filtered.sorted.bam
	""" > ../bsubFiles/filterCapSeqReads.${time}.bsub
done
