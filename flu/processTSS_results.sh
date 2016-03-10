#awk '{ if (NR > 1) {printf("%s\t%d\t%d\t%s\t%.2f\t%.2f\t%.2f\t%.2f\n",$1,$2,$2,$3,($4+1)/($8+1),($5+1)/($9+1),($6+1)/($10+1),($7+1)/($11+1))} else {print "chr\ttss\ttss\tstrand\t6h\t12h\t24h\t48h"}}' /home/si14w/gnearline/flu/txt/tss_counts.txt > /home/si14w/gnearline/flu/txt/tss_counts.foldChange.txt

awk '{ if ($6 == "-") {tss=$3} else if ($6 == "+") {tss=$2} printf("%s\t%d\t%d\t%s\t%s\t%s\n", $1,tss,tss,$4,$5,$6) }' ../txt/hg19_one_tr_per_gene.bed > ../txt/hg19_refseq_genes_plus_snRNA_snoRNA.tss.bed
awk '{ if ($6 == "-") { tss=$3} else if ($6 == "+") { tss=$2} printf("%s\t%d\t%d\t%s_%s_%d_%d\t%s\t%s\n",$1,tss,tss,$4,$1,$2,$3,$5,$6)}' ../txt/hg19_snRNA_and_snoRNA.tmp.bed >> ../txt/hg19_refseq_genes_plus_snRNA_snoRNA.tss.bed
sort -k1,1 -k2,2n ../txt/hg19_refseq_genes_plus_snRNA_snoRNA.tss.bed > ../txt/hg19_refseq_genes_plus_snRNA_snoRNA.tss.sorted.bed

for time in {"A","B","C","D","E","F","G","H"}
do
	echo $time
	awk '{ if ($1 != "track")
		{
			if ($4 > 0)
				printf("%s\t%d\t%d\t%s@%d@%s\t%.2f\t+\n",$1,$2,$3,$1,$2,$3,$4)
			else if ($4 < 0)
				printf("%s\t%d\t%d\t%s@%d@%s\t%.2f\t-\n",$1,$2,$3,$1,$2,$3,(-1*$4))
		}
	}' /home/si14w/gnearline/flu/txt/${time}_normalized.bedGraph | sort -k1,1 -k2,2n > /home/si14w/gnearline/flu/txt/${time}_normalized.bed
	top25_thresh=`awk '{ print $5}' /home/si14w/gnearline/flu/txt/${time}_normalized.bed | Rscript $HOME/gnearline/hdc/scripts/getQuantile.R | cut -d" " -f4`
	awk -vthresh=$top25_thresh '{ if ($5 > thresh) print }' /home/si14w/gnearline/flu/txt/${time}_normalized.bed >  /home/si14w/gnearline/flu/txt/${time}_normalized.tpm_above_75thperc.bed
	bedtools merge -i  /home/si14w/gnearline/flu/txt/${time}_normalized.tpm_above_75thperc.bed -s -d 20 -c 4,5 -o collapse,collapse | awk -f pickTss.awk | sort -k1,1 -k2,2n > /home/si14w/gnearline/flu/txt/${time}_normalized.TSS.sorted.bed
# Use below line if we need to account for closeness to annotated TSS. Right now, only picking TSS with most tag counts.
#	awk -vthresh=$top25_thresh '{ if ($5 > thresh) print }' /home/si14w/gnearline/flu/txt/${time}_normalized.bed | bedtools closest -a stdin -b ~/gnearline/common_scripts/hg19_refseq_rsem_corrected_tss.sorted.bed -d | bedtools merge -i stdin -s -d 100 -c 4,5,13,10 -o collapse,collapse,collapse,distinct 
done

