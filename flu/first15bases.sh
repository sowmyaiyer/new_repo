awk '{ if ($6 == "-") {start=$3-15;end=$3} else if ($6 == "+") {start=$2;end=$2+15} printf("%s\t%d\t%d\t%s\t%s\t%s\n", $1,start,end,$4,$5,$6) }' ../txt/hg19_ucsc_genes_and_flu.bed > ../txt/hg19_ucsc_genes_and_flu.5prime_plus15.allRNA.bed
#sort -k4,4 hg19_snRNA.bed > hg19_snRNA.sorted_by_genename.bed
#sort -k1,1 hg19_ens_to_genename.bed > hg19_ens_to_genename.sorted_by_ens_genename.bed
#join hg19_snRNA.sorted_by_genename.bed hg19_ens_to_genename.sorted_by_ens_genename.bed -1 4 -2 1| awk -F" " '{ print $2"\t"$3"\t"$4"\t"$NF"\t"0"\t"$6}' | sort -k1,1 -k2,2n > ../txt/hg19_snRNA.sorted.bed
awk '{ if ($6 == "-") {start=$3-15;end=$3} else if ($6 == "+") {start=$2;end=$2+15} printf("%s\t%d\t%d\t%s\t%s\t%s\n", $1,start,end,$4,$5,$6) }' ../txt/hg19_snRNA.sorted.bed >> ../txt/hg19_ucsc_genes_and_flu.5prime_plus15.allRNA.bed
