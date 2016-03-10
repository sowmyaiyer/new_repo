## snRNA and snoRNA have one exon each, not spliced. Checked that txStart and txEnd is the same as exon start and exon end.
#grep -v '^#' hg19_snRNA_and_snoRNA | awk '{ printf("%s\t%d\t%d\t%s_%s_%s\t1\t%s\n",$3,$5,$6,$17,$19,$13,$4) }' > hg19_snRNA_and_snoRNA.exons.bed
#cat /project/umw_garberlab/edonnard/dc_project/intron_analysis/reference/exons.bed hg19_snRNA_and_snoRNA.exons.bed > hg19_exons_elisa_with_snRNA_and_snoRNA.bed
#cp /project/umw_garberlab/edonnard/dc_project/intron_analysis/reference/introns.bed ./hg19_introns_elisa.bed
echo "sample tss_count_in_exon tss_count_in_intron tss_count_in_intergenic_regions total_tss_count" | sed 's/ /\t/g' >  ../txt/TSS_annotation.txt
for time in {"A","B","C","D","E","F","G","H"}
do
	tss_count_in_exon=`bedtools intersect -s -a ../txt/${time}_normalized.tpm_above_75thperc.bed -b ../txt/hg19_exons_elisa_with_snRNA_and_snoRNA.bed | wc -l`
	tss_count_in_intron=`bedtools intersect -s -a ../txt/${time}_normalized.tpm_above_75thperc.bed -b ../txt/hg19_introns_elisa.bed | wc -l`
	tss_count_in_intergenic_regions=`bedtools intersect -a ../txt/${time}_normalized.tpm_above_75thperc.bed -b ../txt/hg19_intergenic.elisa.bed | wc -l`
	total_tss_count=`wc -l ../txt/${time}_normalized.tpm_above_75thperc.bed | cut -d" " -f1`
	echo ${time} $tss_count_in_exon $tss_count_in_intron $tss_count_in_intergenic_regions $total_tss_count | sed 's/ /\t/g' >> ../txt/TSS_annotation.txt
	bedtools intersect -s -a ../txt/${time}_normalized.tpm_above_75thperc.bed -b ../txt/hg19_exons_elisa_with_snRNA_and_snoRNA.bed | cut -f5 > ../txt/tss_tpm_exons_${time}.txt
	bedtools intersect -s -a ../txt/${time}_normalized.tpm_above_75thperc.bed -b ../txt/hg19_introns_elisa.bed | cut -f5  > ../txt/tss_tpm_introns_${time}.txt
	bedtools intersect -a ../txt/${time}_normalized.tpm_above_75thperc.bed -b ../txt/hg19_intergenic.elisa.bed | cut -f5 > ../txt/tss_tpm_intergenic_${time}.txt
done
