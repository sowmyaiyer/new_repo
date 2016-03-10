for time in {"A","B","C","D","E","F","G","H"}
do
	echo $time
	sort -k4,4  /home/si14w/gnearline/flu/txt/112213flucap_${time}.readCounts_150bp.allRNA.genes.bed | awk '{ print $4"\t"$NF}' > /home/si14w/gnearline/flu/txt/112213flucap_${time}.readCountsOnly.genes.txt
	join /home/si14w/gnearline/flu/txt/112213flucap_${time}.readCountsOnly.genes.txt /home/si14w/gnearline/flu/txt/gene_mappings_unique_reads.${time}.txt -1 1 -2 1 >
done
