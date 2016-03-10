for bamFile in /project/umw_garberlab/human/DC/ATAC-Seq/ATAC_ChIP_12_22_15/pipe_E70_E72_merged/alignments/*bam
do
	time=`basename $bamFile | cut -d"." -f1`
	echo """
	
	macs2 callpeak --broad -t $bamFile -f BAM -g hs -n ATAC_broad_${time} --outdir ../HumanDC_ATAC_broad/""" > ../bsubFiles/macs2_broad_nucleoatac_dc_${time}.bsub
done
