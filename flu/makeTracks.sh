for time in {"A","B","C","D","E","F","G","H"}
do
	echo """
	module load IGVTools/2.3.31
	sed '/chr[1-9]_[a-z]/d' ../txt/genes_t_df_final${time}.txt | sed '/chr[1-9][0-9]_[a-z]/d' |  sed '/chrUn/d' | sed 's/_ENSG/-ENSG/g' | awk '{
		if (NR > 1) {
		split(\$1,arr,\"_\")
		split(arr[3],coords,\"@\")
		printf(\"%s\\t%d\\t%d\\t%.2f\\n\", coords[1],coords[2],coords[3],\$2)
		}
		
	}' | bedtools sort -i stdin > ../txt/snatchCounts.${time}.wig
	igvtools toTDF -z 5 ../txt/snatchCounts.${time}.wig ../txt/snatchCounts.${time}.tdf ../txt/hg19_plus_flu.chrom.sizes


	sed '/chr[1-9]_[a-z]/d' ../txt/genes_t_df_final${time}.txt | sed '/chr[1-9][0-9]_[a-z]/d' |  sed '/chrUn/d' | sed 's/_ENSG/-ENSG/g' | awk '{
                if (NR > 1) {
                split(\$1,arr,\"_\")
                split(arr[3],coords,\"@\")
                printf(\"%s\\t%d\\t%d\\t%.2f\\n\", coords[1],coords[2],coords[3],\$3)
                }
        }' | bedtools sort -i stdin >  ../txt/capCounts.${time}.wig
        igvtools toTDF -z 5 ../txt/capCounts.${time}.wig ../txt/capCounts.${time}.tdf ../txt/hg19_plus_flu.chrom.sizes
	""" > ../bsubFiles/makeTracks.${time}.bsub

done
