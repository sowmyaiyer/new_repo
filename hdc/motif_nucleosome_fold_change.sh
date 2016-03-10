while read line
do
	motif_id=`echo $line | awk '{ print $1}'`
	motif_names=`echo $line | awk '{ print $2}'`
	echo $motif_names
	if [[ -f /project/umw_garberlab/edonnard/motifs/human/bedfiles/per_motif/matches_${motif_id}.bed ]]; then
	motif_length=`head -1 /project/umw_garberlab/edonnard/motifs/human/bedfiles/per_motif/matches_$motif_id.bed | awk '{ print ($3-$2)}'`
	echo """
		bedtools slop -i  /project/umw_garberlab/edonnard/motifs/human/bedfiles/per_motif/matches_${motif_id}.bed -b 100 -g ~/gnearline/common_scripts/hg19.sorted.genome | awk 'BEGIN{OFS=\"\\t\"} { print \$1,\$2,\$3,NR}' | bigWigAverageOverBed  /home/si14w/gnearline/hdc/nucleoatac_out/nucleoatac_out_E70_E72_0h.occ.bw stdin ~/gnearline/RECYCLE_BIN/out_${motif_id}_0h.tab -bedOut=/farline/umw_flustore_schiffer/hdc/nucleoatac_out_E70_E72_0h.${motif_id}.bed
		bedtools slop -i  /project/umw_garberlab/edonnard/motifs/human/bedfiles/per_motif/matches_${motif_id}.bed -b 100 -g ~/gnearline/common_scripts/hg19.sorted.genome | awk 'BEGIN{OFS=\"\\t\"} { print \$1,\$2,\$3,NR}' | bigWigAverageOverBed  /home/si14w/gnearline/hdc/nucleoatac_out/nucleoatac_out_E70_E72_2h_LPS.occ.bw stdin ~/gnearline/RECYCLE_BIN/out_${motif_id}_2h_LPS.tab -bedOut=/farline/umw_flustore_schiffer/hdc/nucleoatac_out_E70_E72_2h_LPS.${motif_id}.bed
		sort -k1,1 -k2,2n /farline/umw_flustore_schiffer/hdc/nucleoatac_out_E70_E72_0h.${motif_id}.bed | awk '{ print \$NF}' > /farline/umw_flustore_schiffer/hdc/nucleoatac_out_E70_E72_0h.${motif_id}.mean_scores.txt
		sort -k1,1 -k2,2n /farline/umw_flustore_schiffer/hdc/nucleoatac_out_E70_E72_2h_LPS.${motif_id}.bed | awk '{ print \$NF}' > /farline/umw_flustore_schiffer/hdc/nucleoatac_out_E70_E72_2h_LPS.${motif_id}.mean_scores.txt
		rm /farline/umw_flustore_schiffer/hdc/nucleoatac_out_E70_E72_0h.${motif_id}.bed /farline/umw_flustore_schiffer/hdc/nucleoatac_out_E70_E72_2h_LPS.${motif_id}.bed
""" > ../bsubFiles/mean_nucleosome_around_$motif_id.bsub
	fi
done</project/umw_garberlab/edonnard/motifs/human/motif_to_TF_map.tsv
