#cat /project/umw_garberlab/human/DC/ChIP-Seq/ChIP_seq_08_05_2015/pipeline/seqmapping/chip/E01_6d_0h_H3K4me3.1.fastq /project/umw_garberlab/human/DC/ChIP-Seq/ChIP_seq_08_05_2015/pipeline/seqmapping/chip/E01_6d_0h_H3K4me3.2.fastq > ~/gnearline/RECYCLE_BIN/E01_6d_0h_H3K4me3.se.fastq
#bwa index -p ../txt/hg19_repeat_masked.fa ../txt/hg19.fa.align.gz
#bwa aln /share/data/umw_biocore/Genomes/human/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.6.0/genome.fa ~/gnearline/RECYCLE_BIN/E01_6d_0h_H3K4me3.se.fastq > ../bwa_out/E01_6d_0h_H3K4me3.se.sai
#bwa samse /share/data/umw_biocore/Genomes/human/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.6.0/genome.fa ../bwa_out/E01_6d_0h_H3K4me3.se.sai ~/gnearline/RECYCLE_BIN/E01_6d_0h_H3K4me3.se.fastq > ../bwa_out/E01_6d_0h_H3K4me3.se.sam
bwa aln ../txt/hg19.masked.fa ~/gnearline/RECYCLE_BIN/E01_6d_0h_H3K4me3.se.fastq > ../bwa_out/E01_6d_0h_H3K4me3.se.repeat_masked.sai
bwa samse ../txt/hg19.masked.fa ../bwa_out/E01_6d_0h_H3K4me3.se.repeat_masked.sai ~/gnearline/RECYCLE_BIN/E01_6d_0h_H3K4me3.se.fastq > ../bwa_out/E01_6d_0h_H3K4me3.se.repeat_masked.sam
samtools view -bS -o ../bwa_out/E01_6d_0h_H3K4me3.se.repeat_masked.bam ../bwa_out/E01_6d_0h_H3K4me3.se.repeat_masked.sam
