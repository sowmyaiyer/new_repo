awk '{ print $4}' ../txt/hg19_ucsc_genes_and_flu.smallRNAs.bed > ../txt/smallRNAs.txt
while read transcript
do
	echo $transcript
	sed "/${transcript}/d" ../txt/hg19_ucsc_genes_and_flu.gtf > tmp
	grep $transcript tmp
	mv tmp ../txt/hg19_ucsc_genes_and_flu.gtf
done<../txt/smallRNAs.txt
