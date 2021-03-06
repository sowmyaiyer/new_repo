while read line
do
	motif_id=`echo $line | awk '{ print $1}'`
	motif_names=`echo $line | awk '{ print $2}'`
	echo $motif_names
	if [[ -f /project/umw_garberlab/edonnard/motifs/human/bedfiles/per_motif/matches_${motif_id}.bed ]]; then
	motif_length=`head -1 /project/umw_garberlab/edonnard/motifs/human/bedfiles/per_motif/matches_$motif_id.bed | awk '{ print ($3-$2)}'`
	echo """
	awk '{ if (\$6 == \"-\") print \$1\"\\t\"\$2\"\\t\"\$3\"\\tmatch\"NR\"\\t\"\$5\"\\t\"\$6}' /project/umw_garberlab/edonnard/motifs/human/bedfiles/per_motif/matches_$motif_id.bed > /farline/umw_flustore_schiffer/hdc/matches_$motif_id.minus.bed
	awk '{ if (\$6 == \"+\") print \$1\"\\t\"\$2\"\\t\"\$3\"\\tmatch\"NR\"\\t\"\$5\"\\t\"\$6 }' /project/umw_garberlab/edonnard/motifs/human/bedfiles/per_motif/matches_$motif_id.bed >  /farline/umw_flustore_schiffer/hdc/matches_$motif_id.plus.bed
	for time in {\"0h\",\"30min_LPS\",\"2h_LPS\",\"4h_LPS\",\"24h_LPS\"}
	do
		bedtools slop -i  /farline/umw_flustore_schiffer/hdc/matches_$motif_id.minus.bed -b 100 -g ~/gnearline/common_scripts/hg19.sorted.genome | bedtools makewindows -b stdin -w 1 -reverse -i srcwinnum | bigWigAverageOverBed  ../txt/merged_E81_\${time}.cutadapt.sorted.no_dups.filt.nucSignal.centRL74.se.bam.bw stdin ~/gnearline/RECYCLE_BIN/out_${motif_id}_\${time}.minus.tab -bedOut=/farline/umw_flustore_schiffer/hdc/merged_E81_\${time}.${motif_id}.nuclSignal.minus.bed
		bedtools slop -i  /farline/umw_flustore_schiffer/hdc/matches_$motif_id.minus.bed -b 100 -g ~/gnearline/common_scripts/hg19.sorted.genome | bedtools makewindows -b stdin -w 1 -reverse -i srcwinnum | bigWigAverageOverBed  ../txt/merged_E81_\${time}.cutadapt.sorted.no_dups.filt.nucFree.centRL40.se.bam.bw stdin ~/gnearline/RECYCLE_BIN/out_${motif_id}_\${time}.minus.tab -bedOut=/farline/umw_flustore_schiffer/hdc/merged_E81_\${time}.${motif_id}.nucFree.minus.bed
		rm ~/gnearline/RECYCLE_BIN/out_${motif_id}_\${time}.minus.tab
		bedtools slop -i  /farline/umw_flustore_schiffer/hdc/matches_$motif_id.plus.bed -b 100 -g ~/gnearline/common_scripts/hg19.sorted.genome | bedtools makewindows -b stdin -w 1 -i srcwinnum | bigWigAverageOverBed  ../txt/merged_E81_\${time}.cutadapt.sorted.no_dups.filt.nucSignal.centRL74.se.bam.bw stdin ~/gnearline/RECYCLE_BIN/out_${motif_id}_\${time}.plus.tab -bedOut=/farline/umw_flustore_schiffer/hdc/merged_E81_\${time}.${motif_id}.nuclSignal.plus.bed
		bedtools slop -i  /farline/umw_flustore_schiffer/hdc/matches_$motif_id.plus.bed -b 100 -g ~/gnearline/common_scripts/hg19.sorted.genome | bedtools makewindows -b stdin -w 1 -i srcwinnum | bigWigAverageOverBed  ../txt/merged_E81_\${time}.cutadapt.sorted.no_dups.filt.nucFree.centRL40.se.bam.bw stdin ~/gnearline/RECYCLE_BIN/out_${motif_id}_\${time}.plus.tab -bedOut=/farline/umw_flustore_schiffer/hdc/merged_E81_\${time}.${motif_id}.nucFree.plus.bed
		rm ~/gnearline/RECYCLE_BIN/out_${motif_id}_\${time}.plus.tab
                sed 's/_/\t/g'  /farline/umw_flustore_schiffer/hdc/merged_E81_\${time}.${motif_id}.nuclSignal.minus.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' > /farline/umw_flustore_schiffer/hdc/merged_E70_E72.\${time}.${motif_id}.nuclSignal.all_signals.txt
                sed 's/_/\t/g'  /farline/umw_flustore_schiffer/hdc/merged_E81_\${time}.${motif_id}.nuclSignal.plus.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' >> /farline/umw_flustore_schiffer/hdc/merged_E70_E72.\${time}.${motif_id}.nuclSignal.all_signals.txt
                sed 's/_/\t/g'  /farline/umw_flustore_schiffer/hdc/merged_E81_\${time}.${motif_id}.nucFree.minus.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' > /farline/umw_flustore_schiffer/hdc/merged_E70_E72.\${time}.${motif_id}.nucFree.all_signals.txt
                sed 's/_/\t/g'  /farline/umw_flustore_schiffer/hdc/merged_E81_\${time}.${motif_id}.nucFree.plus.bed | sort -k4,4 -k5,5n | bedtools groupby -g 4 -c 6 -o collapse | sed 's/,/\t/g' >> /farline/umw_flustore_schiffer/hdc/merged_E70_E72.\${time}.${motif_id}.nucFree.all_signals.txt
		rm  /farline/umw_flustore_schiffer/hdc/merged_E81_\${time}.${motif_id}.*.*us.bed
	done
	rm  /farline/umw_flustore_schiffer/hdc/matches_$motif_id.minus.bed  /farline/umw_flustore_schiffer/hdc/matches_$motif_id.plus.bed 
""" > ../bsubFiles/nucleosome_around_$motif_id.bsub
	fi
done</project/umw_garberlab/edonnard/motifs/human/motif_to_TF_map.tsv
