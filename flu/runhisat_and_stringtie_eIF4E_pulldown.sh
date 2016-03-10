fastq_dir="/farline/umw_flustore_schiffer/eIF4E_pulldown_RNAseq/fastq_untrimmed"
output_dir="/project/umw_flustore_schiffer/fluseq/eIF4E_pulldown_RNAseq/hisat2_out"
while read sample
do
	echo """
	module load cutadapt/1.9
	module load IGVTools/2.3.31
	module load stringtie/1.1.2
	cutadapt -u -75 -o ${fastq_dir}/${sample}.R1.trimmed.fastq.gz ${fastq_dir}/${sample}.R1.fastq.gz
	hisat2 -p 12 --rna-strandness R -x /project/umw_flustore_schiffer/fluseq/hg19_and_fluWithSnp_combined_hisat_index -U ${fastq_dir}/${sample}.R1.trimmed.fastq.gz -S ${output_dir}/hisat2_out_${sample}.sam
	samtools view -bS -F 4 -o ${output_dir}/hisat2_out_${sample}.bam  ${output_dir}/hisat2_out_${sample}.sam
        samtools sort -T ${output_dir}/hisat2_out_${sample}.sorted -o ${output_dir}/hisat2_out_${sample}.sorted.bam  ${output_dir}/hisat2_out_${sample}.bam
        samtools index ${output_dir}/hisat2_out_${sample}.sorted.bam
	samtools view -b -F 16 -o ${output_dir}/hisat2_out_${sample}.sorted.flu_vRNA.bam ${output_dir}/hisat2_out_${sample}.sorted.bam 
	samtools view -b -f 16 -o ${output_dir}/hisat2_out_${sample}.sorted.flu_mRNA.bam ${output_dir}/hisat2_out_${sample}.sorted.bam 
	rm   ${output_dir}/hisat2_out_${sample}.bam ${output_dir}/hisat2_out_${sample}.sam
	mkdir -p ${output_dir}/tdf
	igvtools count -z 5 -w 10 -e 200 ${output_dir}/hisat2_out_${sample}.sorted.bam  ${output_dir}/tdf/hisat2_out_${sample}.sorted.tdf ../txt/hg19_plus_flu.chrom.sizes
	igvtools count -z 5 -w 10 -e 200 ${output_dir}/hisat2_out_${sample}.sorted.flu_mRNA.bam  ${output_dir}/tdf/hisat2_out_${sample}.sorted.flu_mRNA.tdf ../txt/hg19_plus_flu.chrom.sizes
	igvtools count -z 5 -w 10 -e 200 ${output_dir}/hisat2_out_${sample}.sorted.flu_vRNA.bam  ${output_dir}/tdf/hisat2_out_${sample}.sorted.flu_vRNA.tdf ../txt/hg19_plus_flu.chrom.sizes
	mkdir -p ${output_dir}/stringtie_out
	stringtie ${output_dir}/hisat2_out_${sample}.sorted.bam -p 12 -G $HOME/gnearline/flu/txt/biocore_ucsc_plus_flu.gtf -o ${output_dir}/stringtie_out/stringtie_out_${sample}.gtf
	""" > ../bsubFiles/hisat2_${sample}.bsub
done<samples.txt
