samtools flagstat ../txt/GM12878_H3K4me3_Bernstein.bam > ../txt/GM12878_H3K4me3_Bernstein.flagstat
mappedReads_se=`awk '{ if (NR ==5) print $1}' ../txt/GM12878_H3K4me3_Bernstein.flagstat`
scalingFactor_se=`echo $mappedReads_se | awk -vmappedReads=$mappedReads_se 'BEGIN{print 1000000/mappedReads}'`
samtools sort -T ../txt/GM12878_H3K4me3_Bernstein.sorted -o ../txt/GM12878_H3K4me3_Bernstein.sorted.bam ../txt/GM12878_H3K4me3_Bernstein.bam
bedtools genomecov -ibam ../txt/GM12878_H3K4me3_Bernstein.sorted.bam -g ~/gnearline/common_scripts/hg19.genome -bga -scale ${scalingFactor_se} > ../bwa_out/GM12878_H3K4me3_Bernstein.se.bedGraph
bedGraphToBigWig ../bwa_out/GM12878_H3K4me3_Bernstein.se.bedGraph ~/gnearline/common_scripts/hg19.genome ../bwa_out/GM12878_H3K4me3_Bernstein.se.bw
