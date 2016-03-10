for time in {"A","B","C","D","E","F","G","H"}
do
	echo $time
#	mv ../bowtie_out/112213flucap_${time}.sorted.bam ../bowtie_out/112213flucap_${time}.softclip.sorted.bam
	echo """
	samtools view ../bowtie_out/112213flucap_${time}.softclip.sorted.bam | grep "XM:i:0" | awk '{ print \$1}' | sort | uniq -c | awk '{ print \$1}' >  ../txt/histogram_multireads_${time}.txt
#	awk '{ print \$1}' ../txt/multiread_dist_capseq.${time}.txt > ../txt/histogram_multireads_${time}.txt
	""" > ../bsubFiles/multiread_numbers_${time}.bsub
done
