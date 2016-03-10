for time in {"A","B","C","D","E","F","G","H"}
do
	echo $time
	samtools view -F 4 ../bowtie_out/112213flucap_rRNA_${time}.sorted.bam | awk '{ print $3}' | sort | uniq -c > ../txt/rRNA_reads_time_${time}.txt
	total_rRNA_reads=`awk 'BEGIN{sum=0}{sum=sum+$1}END{print sum}' ../txt/rRNA_reads_time_${time}.txt`
	echo "total rRNA reads = ${total_rRNA_reads}" >> ../txt/rRNA_reads_time_${time}.txt
	totalReads=`samtools view -c ../bowtie_out/112213flucap_${time}.sorted.bam`
	echo "total reads = ${totalReads}" >>  ../txt/rRNA_reads_time_${time}.txt
done
