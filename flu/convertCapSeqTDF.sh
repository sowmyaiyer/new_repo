for time in {"A","B","C","D","E","F","G","H"}
do
		echo """
		module load IGVTools/2.3.31
		igvtools count -z 5 -w 25 -e 100 ../bowtie_out/112213flucap_${time}.sorted.bam ../bowtie_out/112213flucap_${time}.sorted.tdf ../txt/hg19_plus_flu.chrom.sizes""" > ../bsubFiles/capseq_tdf_$time.bsub
done
