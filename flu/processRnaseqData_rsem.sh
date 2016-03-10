#awk '{ if (($3 -$2) < 200) print }'  ../txt/hg19_ucsc_genes_and_flu.bed > ../txt/hg19_ucsc_genes_and_flu.smallRNAs.bed


#grep -vwF -f ../txt/smallRNAs.txt ../txt/hg19_ucsc_genes_and_flu.gtf > ../txt/hg19_ucsc_genes_and_flu.longer_than_200.gtf
# Had to manually remove a blank line in ../txt/hg19_ucsc_genes_and_flu.longer_than_200.gtf. Why?

for time in {"6h","12h","24h","48h","mock12h"}
do
echo """#BSUB -w done(4699985)
module unload R/3.2.0
module load RSEM/1.2.11
#rsem-prepare-reference --gtf ../txt/hg19_ucsc_genes_and_flu.longer_than_200.gtf --transcript-to-gene-map ../txt/knownIsoforms_hg19_and_flu.txt --bowtie2 --bowtie2-path /share/pkg/bowtie2/2-2.1.0 ../hg19_and_fluWithSnp_combined.fa ../rsem_reference_files/hg19_and_fluWithSnp_combined 
rsem-calculate-expression --bowtie2 --bowtie2-path /share/pkg/bowtie2/2-2.1.0 -p 8 --forward-prob 0 --fragment-length-mean 180 --fragment-length-sd 50 /farline/umw_flustore_schiffer/rnaseq/04SEP15.A549_B59_${time}.fastq.gz ../rsem_reference_files/hg19_and_fluWithSnp_combined ../txt/04SEP15.A549_B59_totalRNA_${time}
""" > ../bsubFiles/rsem_flu_${time}.bsub
done
