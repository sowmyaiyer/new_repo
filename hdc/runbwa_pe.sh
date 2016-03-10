bwa aln ../txt/hg19.masked.fa ../txt/E01_6d_0h_H3K4me3.1.fastq > ../bwa_out/E01_6d_0h_H3K4me3.1.masked.sai
bwa aln ../txt/hg19.masked.fa ../txt/E01_6d_0h_H3K4me3.2.fastq > ../bwa_out/E01_6d_0h_H3K4me3.2.masked.sai
bwa sampe ../txt/hg19.masked.fa ../bwa_out/E01_6d_0h_H3K4me3.1.masked.sai ../bwa_out/E01_6d_0h_H3K4me3.2.masked.sai ../txt/E01_6d_0h_H3K4me3.1.fastq ../txt/E01_6d_0h_H3K4me3.2.fastq > ../bwa_out/E01_6d_0h_H3K4me3.masked.pe.sam
samtools view -bS -o ../bwa_out/E01_6d_0h_H3K4me3.masked.pe.bam ../bwa_out/E01_6d_0h_H3K4me3.masked.pe.sam
