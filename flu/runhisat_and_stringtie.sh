for time in {"6h","12h","24h","48h","mock12h"}
do
	echo """
	module load stringtie/1.1.2
	module load IGVTools/2.3.31
	hisat2 -p 12 --rna-strandness R -x /project/umw_flustore_schiffer/fluseq/hg19_and_fluWithSnp_combined_hisat_index -U /farline/umw_flustore_schiffer/rnaseq/04SEP15.A549_B59_${time}_all.fastq.gz -S /project/umw_flustore_schiffer/fluseq/old_rnaseq/hisat2_out_${time}.sam
	samtools view -bS -F 4 -o /project/umw_flustore_schiffer/fluseq/old_rnaseq/hisat2_out_${time}.bam  /project/umw_flustore_schiffer/fluseq/old_rnaseq/hisat2_out_${time}.sam
        samtools sort -T /project/umw_flustore_schiffer/fluseq/old_rnaseq/hisat2_out_${time}.sorted -o /project/umw_flustore_schiffer/fluseq/old_rnaseq/hisat2_out_${time}.sorted.bam  /project/umw_flustore_schiffer/fluseq/old_rnaseq/hisat2_out_${time}.bam
        samtools index /project/umw_flustore_schiffer/fluseq/old_rnaseq/hisat2_out_${time}.sorted.bam
	rm  /project/umw_flustore_schiffer/fluseq/old_rnaseq/hisat2_out_${time}.sam  /project/umw_flustore_schiffer/fluseq/old_rnaseq/hisat2_out_${time}.bam
	mkdir  /project/umw_flustore_schiffer/fluseq/old_rnaseq/tdf
	igvtools count -z 5 -w 10 -e 200 /project/umw_flustore_schiffer/fluseq/old_rnaseq/hisat2_out_${time}.sorted.bam  /project/umw_flustore_schiffer/fluseq/old_rnaseq/tdf/hisat2_out_${time}.sorted.tdf ../txt/hg19_plus_flu.chrom.sizes
	mkdir /project/umw_flustore_schiffer/fluseq/old_rnaseq/stringtie_out
	stringtie  /project/umw_flustore_schiffer/fluseq/old_rnaseq/hisat2_out_${time}.sorted.bam -p 12 -G $HOME/gnearline/flu/txt/biocore_ucsc_plus_flu.gtf -o  /project/umw_flustore_schiffer/fluseq/old_rnaseq/stringtie_out/stringtie_out_${time}.gtf
	""" > ../bsubFiles/hisat2_${time}.bsub
done
