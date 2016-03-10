# ./processChipNew_pe_JK_vs_DC_preprocess.sh
for sample in  `awk '{ print $NF}' ../txt/chip01132016_sampleInfo.txt`
do
	echo $sample
	echo """
	module load IGVTools/2.3.31
	igvtools count -z 5 -w 25 -e 250 ../bwa_out/${sample}.pe.properly_paired.sorted.bam  ../bwa_out/${sample}.pe.properly_paired.tdf  hg19
	""" > ../bsubFiles/new_JK_vs_DC_generateTDF_${sample}.bsub
done
