# Get reads whose fragment sizes are > 146 nucleosomes. Make read-depth-normalized signal file
for time in {"0h","30min_LPS","2h_LPS","4h_LPS","24h_LPS"}
do
	echo $time
	nuclSignal_bam=/project/umw_garberlab/human/DC/ATAC-Seq/ATAC_E81_merge/pipe2/tssPileups/Centered/merged_E81_${time}.cutadapt.sorted.no_dups.filt.nucSignal.centRL74.se.bam
	nuclFreeSignal_bam=/project/umw_garberlab/human/DC/ATAC-Seq/ATAC_E81_merge/pipe2/tssPileups/Centered/merged_E81_${time}.cutadapt.sorted.no_dups.filt.nucFree.centRL40.se.bam
	totalReads_nuclSignal=`samtools view -c ${nuclSignal_bam}`
	totalReads_nuclFree=`samtools view -c ${nuclFreeSignal_bam}`
	scalingFactor_nuclSignal=`echo $totalReads_nuclSignal | awk -vmappedReads=$totalReads_nuclSignal 'BEGIN{print 1000000/mappedReads}'`
	scalingFactor_nuclFree=`echo $totalReads_nuclFree | awk -vmappedReads=$totalReads_nuclFree 'BEGIN{print 1000000/mappedReads}'`
        bedtools genomecov -ibam ${nuclSignal_bam} -g ~/gnearline/common_scripts/hg19.sorted.genome -bga -scale ${scalingFactor_nuclSignal} > ../txt/merged_E81_${time}.cutadapt.sorted.no_dups.filt.nucSignal.centRL74.se.bam.bedGraph
        bedtools genomecov -ibam ${nuclFreeSignal_bam} -g ~/gnearline/common_scripts/hg19.sorted.genome -bga -scale ${scalingFactor_nuclFree} > ../txt/merged_E81_${time}.cutadapt.sorted.no_dups.filt.nucFree.centRL40.se.bam.bedGraph
        bedGraphToBigWig  ../txt/merged_E81_${time}.cutadapt.sorted.no_dups.filt.nucSignal.centRL74.se.bam.bedGraph ~/gnearline/common_scripts/hg19.sorted.genome  ../txt/merged_E81_${time}.cutadapt.sorted.no_dups.filt.nucSignal.centRL74.se.bam.bw
        bedGraphToBigWig  ../txt/merged_E81_${time}.cutadapt.sorted.no_dups.filt.nucFree.centRL40.se.bam.bedGraph ~/gnearline/common_scripts/hg19.sorted.genome  ../txt/merged_E81_${time}.cutadapt.sorted.no_dups.filt.nucFree.centRL40.se.bam.bw
	rm ../txt/merged_E81_${time}.cutadapt.sorted.no_dups.filt.nucSignal.centRL74.se.bam.bedGraph ../txt/merged_E81_${time}.cutadapt.sorted.no_dups.filt.nucFree.centRL40.se.bam.bedGraph
done
