while read line
do
	motif_code=`echo $line | awk '{ print $1}'`
	motif_tf=`echo $line | awk '{ print $2}'`
	echo $motif_code
	bedtools closest -a /project/umw_garberlab/edonnard/motifs/human/bedfiles/per_motif/matches_${motif_code}.bed -b /home/si14w/gnearline/hdc/nucleoatac_out/nucleoatac_out_E70_E72_0h.nucmap_combined.bed.gz -d | awk '{ print $NF}' > ../txt/${motif_code}.E70_E72_0h.txt
done</project/umw_garberlab/edonnard/motifs/human/motif_to_TF_map.tsv
