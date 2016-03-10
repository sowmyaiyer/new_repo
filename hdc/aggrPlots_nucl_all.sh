for motif_info_file in /home/si14w/gnearline/hdc/txt/motifs_*
do
	#echo """Rscript aggrPlots_nucl_all.R ${motif_info_file}""" > ../bsubFiles/aggrPlots_nucl_all_`basename ${motif_info_file}`.bsub
	echo """Rscript aggrPlots_nucl_all.raw.R ${motif_info_file}""" > ../bsubFiles/aggrPlots_nucl_all_raw_`basename ${motif_info_file}`.bsub
done
