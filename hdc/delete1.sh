for sample in  `awk '{ print $1}' ../txt/chipnew_sampleInfo.txt | sort | uniq`
do
        echo $sample

	mappedReads_se=`awk '{ if (NR ==5) print $1}' ../bwa_out/${sample}.se.flagstat`
        mappedReads_pe=`awk '{ if (NR ==5) print $1}' ../bwa_out/${sample}.pe.flagstat`

	mappedReadsInPromoter_pe=`samtools view -F 0x4 -b ../bwa_out/${sample}.pe.sorted.bam | bedtools bamtobed -i stdin | bedtools intersect -a stdin -b ../txt/hg19_promoter.bed -u | wc -l`
        mappedReadsInPromoter_se=`samtools view -F 0x4 -b ../bwa_out/${sample}.se.sorted.bam | bedtools bamtobed -i stdin | bedtools intersect -a stdin -b ../txt/hg19_promoter.bed -u | wc -l`

        echo mapped reads in promoter = $mappedReadsInPromoter_pe >> ../bwa_out/${sample}.pe.flagstat
        echo mapped reads in promoter = $mappedReadsInPromoter_se >> ../bwa_out/${sample}.se.flagstat
done
