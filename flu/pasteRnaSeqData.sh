for time in {"6h","12h","24h","48h","mock12h"}
do
	echo $time
	sort -k1,1 ../txt/04SEP15.A549_B59_totalRNA_${time}.isoforms.sorted_by_gene.txt | awk '{ print $2}' > /home/si14w/gnearline/flu/txt/04SEP15.A549_B59_totalRNA_${time}.isoforms.sorted_by_gene.tpm_only.txt
done
sort -k1,1 ../txt/04SEP15.A549_B59_totalRNA_${time}.isoforms.sorted_by_gene.txt | awk '{ print $1}' > ../txt/geneList.txt
echo gene 6h 12h 24h 48h mock_12h > ../txt/rnaseq_rpm_longRNA.txt
paste ../txt/geneList.txt /home/si14w/gnearline/flu/txt/04SEP15.A549_B59_totalRNA_6h.isoforms.sorted_by_gene.tpm_only.txt /home/si14w/gnearline/flu/txt/04SEP15.A549_B59_totalRNA_12h.isoforms.sorted_by_gene.tpm_only.txt /home/si14w/gnearline/flu/txt/04SEP15.A549_B59_totalRNA_24h.isoforms.sorted_by_gene.tpm_only.txt /home/si14w/gnearline/flu/txt/04SEP15.A549_B59_totalRNA_48h.isoforms.sorted_by_gene.tpm_only.txt /home/si14w/gnearline/flu/txt/04SEP15.A549_B59_totalRNA_mock12h.isoforms.sorted_by_gene.tpm_only.txt >>  ../txt/rnaseq_rpm_longRNA.txt

join ../txt/rnaseq_rpm_longRNA.txt ../txt/hg19_refseq-2.sorted.bed -1 1 -2 1 -a 1 | sed 's/ /\t/g' | awk '{ print $1"\t"$7"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > ../txt/rnaseq_rpm_longRNA.gene_alias.txt
