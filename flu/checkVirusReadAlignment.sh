for time in {"A","B","C","D","E","F","G","H"}
do
	echo $time 
	totalReads=`samtools view -F 4 -c ../bowtie_out/112213flucap_${time}.sorted.bam`
	totalVirusReads=`samtools view -F 4 ../bowtie_out/112213flucap_${time}.sorted.bam | grep -c "B59FULL"`
	echo $time $totalVirusReads $totalReads 
done
