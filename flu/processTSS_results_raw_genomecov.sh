for time in {"A","B","C","D","E","F","G","H"}
do
	bedtools genomecov -ibam ../bowtie_out/112213flucap_${time}.sorted.bam -g ../txt/hg19_plus_flu.genome -d -5 > ../txt/${time}_capseq_genomecov.txt
	awk '{ if ($3 > 5) printf("%s\t%d\t%d\t",)}' ../txt/${time}_capseq_genomecov.txt > /home/si14w/gnearline/flu/txt/${time}_raw.tags_above_5.bed
done
