#       bowtie2-build ../hg19_and_fluWithSnp_combined.fa ../bowtie_out/hg19_and_fluWithSnp_combined
for time in {"A","B","C","D","E","F","G","H"}
do
        echo """
        module load bowtie2/2-2.1.0
        bowtie2 --local -x hg19_and_fluWithSnp_combined -U /project/umw_garberlab/narayanan/112213flucap/rawData_toFluAndtoHuman/${time}_Nocodes.fastq.gz -S ../bowtie_out/112213flucap_${time}.sam
        samtools view -bS -o ../bowtie_out/112213flucap_${time}.bam ../bowtie_out/112213flucap_${time}.sam
        samtools sort -T ../bowtie_out/112213flucap_${time}.sorted -o ../bowtie_out/112213flucap_${time}.sorted.bam ../bowtie_out/112213flucap_${time}.bam
        samtools index ../bowtie_out/112213flucap_${time}.sorted.bam
        rm  ../bowtie_out/112213flucap_${time}.sam  ../bowtie_out/112213flucap_${time}.bam

        bedtools genomecov -ibam ../bowtie_out/112213flucap_${time}.sorted.bam -g ../txt/hg19_plus_flu.genome -d -5 > ../txt/${time}_capseq_genomecov.txt
        awk '{ if (\$3 > 5) printf(\"%s\\t%d\\t%d\\t\",)}' ../txt/${time}_capseq_genomecov.txt > /home/si14w/gnearline/flu/txt/${time}_raw.tags_above_5.bed


	samtools view ../bowtie_out/112213flucap_${time}.sorted.bam | grep "B59FULL" | awk '
        {
                if (\$4 == 1)
                {
                        split(\$6,arr,\"S\");
                        if(arr[1] >= 9 && arr[1] <= 15)
                        printf(\">%s\\n%s\\n\",\$1,substr(\$10,1,arr[1]))
                }
        }' > ../bowtie_out/112213flucap_${time}.flu_softclipped.fa


        ## First annotate each TSS according to which gene, intron vs exon etc or if TSS in intergenic region, closest gene, distance to closest gene
        bedtools intersect -s -a ../txt/${time}_raw.tags_above_5.bed -b ../txt/hg19_introns_and_exons.bed.tmp -wb | bedtools groupby -c 10 -o last -full | awk '{ printf(\"%s\\t%d\\t%d\\t%s_%s_%s_%s_0\\t%.2f\\t%s\\n\",\$1,\$2,\$3,\$13,\$10,\$4,\$10,\$5,\$6)}' | bedtools slop -i stdin -r 15 -l 0 -s -g ~/gnearline/common_scripts/hg19.sorted.genome > ../txt/${time}.TSS_plus_15.bed
        bedtools intersect -a ../txt/${time}_raw.tags_above_5.bed -b ../txt/hg19_intergenic.elisa.bed.tmp -wb | bedtools closest -a stdin -b ../txt/hg19_refseq_genes_plus_snRNA_snoRNA.tss.sorted.bed -d | awk '{ printf(\"%s\\t%d\\t%d\\tintergenic_id%s_%s_%s_%d\\t%.2f\\t%s\\n\",\$1,\$2,\$3,NR,\$4,\$14,\$17,\$5,\$6)}' | bedtools slop -i stdin -r 15 -l 0 -s -g ~/gnearline/common_scripts/hg19.sorted.genome >> ../txt/${time}.TSS_plus_15.bed
        bedtools getfasta -name -fi /share/data/umw_biocore/genome_data/human/hg19full/hg19full.fa -bed  ../txt/${time}.TSS_plus_15.bed -fo ../txt/${time}.TSS_plus_15.fa
        bowtie2-build ../txt/${time}.TSS_plus_15.fa ../txt/${time}.TSS_plus_15.fa
        bowtie2 -p 8 -x ../txt/${time}.TSS_plus_15.fa -N 0 --norc -U ../bowtie_out/112213flucap_${time}.flu_softclipped.fa --all -f -S ../bowtie_out/${time}_capseq_TSS.sam
        samtools view -bS -F 4 ../bowtie_out/${time}_capseq_TSS.sam  | samtools sort -n -T ../bowtie_out/${time}_capseq_TSS.sorted -o ../bowtie_out/${time}_capseq_TSS.sorted.bam -
        samtools view  ../bowtie_out/${time}_capseq_TSS.sorted.bam | awk '{ if (\$4 == 1) print }' > ../bowtie_out/${time}_capseq_TSS.sorted.5prime.sam
        grep "^@" ../bowtie_out/${time}_capseq_TSS.sam > ../bowtie_out/${time}_capseq_TSS.samheader
#        rm ../bowtie_out/${time}_capseq_TSS.sam ../bowtie_out/${time}_capseq_TSS.bam


        sort -k4,4 ../txt/${time}.TSS_plus_15.bed > ../txt/${time}.TSS_plus_15.sorted_by_TSS_Id.bed
        awk '{ print \$1}' ../bowtie_out/${time}_capseq_TSS.sorted.5prime.sam | sort | uniq -c | awk '{ if (\$1 == 1) print \$2 }' | sort > ../txt/${time}_capseq_TSS.unique_mapping_reads.txt
        sort -k1,1 ../bowtie_out/${time}_capseq_TSS.sorted.5prime.sam > ../bowtie_out/${time}_capseq_TSS.sorted.5prime.sorted_by_readId.sam
        join  ../bowtie_out/${time}_capseq_TSS.sorted.5prime.sorted_by_readId.sam  ../txt/${time}_capseq_TSS.unique_mapping_reads.txt -1 1 -2 1 | sed 's/ /\t/g' | awk '{ print \$1\"\\t\"\$3}'| sort -k2,2 | bedtools groupby -g 2 -c 1 -o count > ../txt/${time}_capseq_TSS_unique_mappers.incomplete.${time}.txt
        echo "gene capReadCounts snatchReadCounts" | sed 's/ /\t/g' >  ../txt/${time}_capseq_TSS_capCounts_and_unique_snatch_counts.txt
        join ../txt/${time}.TSS_plus_15.sorted_by_TSS_Id.bed ../txt/${time}_capseq_TSS_unique_mappers.incomplete.${time}.txt -1 4 -2 1 -a 1 | sed 's/ /\t/g' | awk '{ if (NF == 6 ) print \$1\"\\t\"\$5\"\\t0\"; else print \$1\"\\t\"\$5\"\\t\"\$7}' >> ../txt/${time}_capseq_TSS_capCounts_and_unique_snatch_counts.txt
        awk '{ print \$1}' ../bowtie_out/${time}_capseq_TSS.sorted.5prime.sam | sort | uniq -c | awk '{ if (\$1 > 1 && \$1 <= 20 ) print \$2 }' | sort > ../txt/${time}_capseq_TSS.multi_mapping_reads.txt
        join  ../bowtie_out/${time}_capseq_TSS.sorted.5prime.sorted_by_readId.sam  ../txt/${time}_capseq_TSS.multi_mapping_reads.txt -1 1 -2 1 | sed 's/ /\t/g' | sort -k1,1 | bedtools groupby -g 1 -c 3,3 -o count,distinct > ../txt/${time}_capseq_TSS.multi_mapping_read_counts.txt
        echo "id multiplicity genenames" | sed 's/ /\t/g' > ../txt/${time}_capseq_TSS.multi_mapping_read_counts.proper.txt
        sort -k3,3 ../txt/${time}_capseq_TSS.multi_mapping_read_counts.txt | bedtools groupby -i stdin -g 3 -c 1 -o count | awk '{ print \"line_\"NR\"\\t\"\$2\"\\t\"\$1}' | sort -k2,2nr >>  ../txt/${time}_capseq_TSS.multi_mapping_read_counts.proper.txt
        Rscript fluSnatch.R /home/si14w/gnearline/flu/txt/${time}_capseq_TSS_capCounts_and_unique_snatch_counts.txt   /home/si14w/gnearline/flu/txt/${time}_capseq_TSS.multi_mapping_read_counts.proper.txt ${time}


	

        """ > ../bsubFiles/flu_capseq_bowtie_${time}.bsub
done

#awk '{ print $0"\tintron"}' hg19_introns_elisa.bed > hg19_introns_elisa.bed.tmp
#awk '{ print $0"\texon"}' hg19_exons_elisa_with_snRNA_and_snoRNA.bed > hg19_exons_elisa_with_snRNA_and_snoRNA.bed.tmp
#awk '{ print $0"\tintergenic"}' hg19_intergenic.elisa.bed > hg19_intergenic.elisa.bed.tmp
#cat hg19_introns_elisa.bed.tmp hg19_exons_elisa_with_snRNA_and_snoRNA.bed.tmp > hg19_introns_and_exons.bed.tmp

for time in {"A","B","C","D"}
do
	echo ${time}
	awk -vtime=${time} '{printf("%s\t%d\t%d\t%s\t%s_%.2f\t%s\n",$1,$2,$3,$4,time,$5,$6) }' ../txt/${time}.TSS_plus_15.bed >> ../txt/TSS_plus_15.combined_time.bed
done


# The following code is to try to map snatches to TSSs across all time points. Only snatches from "C" (24h) were mapped.

bedtools sort -i ../txt/TSS_plus_15.combined_time.bed | bedtools merge -d -16 -s -c 5,4 -o collapse,distinct | awk '{printf("%s\t%d\t%d\t%s|%s\t.\t%s\n",$1,$2,$3,$6,$5,$4) }' >  ../txt/TSS_plus_15.combined_time.annotated.bed
bedtools getfasta -name -fi /share/data/umw_biocore/genome_data/human/hg19full/hg19full.fa -bed ../txt/TSS_plus_15.combined_time.annotated.bed  -fo ../txt/TSS_plus_15.combined_time.fa
module load bowtie2/2-2.1.0
bowtie2 -p 8 -x ../txt/TSS_plus_15.combined_time.fa -N 0 --norc -U ../bowtie_out/112213flucap_C.flu_softclipped.fa --all -f -S ../bowtie_out/C_capseq_TSS.combined_time.sam
samtools view -bS -F 4 ../bowtie_out/C_capseq_TSS.combined_time.sam  | samtools sort -n -T ../bowtie_out/C_capseq_TSS.combined_time.sorted -o ../bowtie_out/C_capseq_TSS.combined_time.sorted.bam -
samtools view  ../bowtie_out/C_capseq_TSS.combined_time.sorted.bam | awk '{ if ($4 == 1) print }' > ../bowtie_out/C_capseq_TSS.sorted.5prime.combined_time.sam
sort -k4,4 ../txt/TSS_plus_15.combined_time.annotated.bed > ../txt/TSS_plus_15.combined_time.annotated.sorted_by_TSS_Id.bed
awk '{ print $1}' ../bowtie_out/C_capseq_TSS.sorted.5prime.combined_time.sam | sort | uniq -c | awk '{ if ($1 == 1) print $2 }' | sort > ../txt/C_capseq_TSS.unique_mapping_reads_combined_time.txt
sort -k1,1 ../bowtie_out/C_capseq_TSS.sorted.5prime.combined_time.sam > ../bowtie_out/C_capseq_TSS.sorted.5prime.combined_time.sorted_by_readId.sam
join  ../bowtie_out/C_capseq_TSS.sorted.5prime.combined_time.sorted_by_readId.sam  ../txt/C_capseq_TSS.unique_mapping_reads_combined_time.txt -1 1 -2 1 | sed 's/ /\t/g' | awk '{ print $1"\t"$3}'| sort -k2,2 | bedtools groupby -g 2 -c 1 -o count > ../txt/C_capseq_TSS_unique_mappers.incomplete.combined_time.txt
echo "gene capReadCounts snatchReadCounts" | sed 's/ /\t/g' >  ../txt/C_capseq_TSS_capCounts_and_unique_snatch_counts_combined_time.txt
join ../txt/TSS_plus_15.combined_time.annotated.sorted_by_TSS_Id.bed ../txt/C_capseq_TSS_unique_mappers.incomplete.combined_time.txt -1 4 -2 1 -a 1 | sed 's/ /\t/g' | awk '{ if (NF == 6 ) print $1"\t"$5"\t0"; else print $1"\t"$5"\t"$7}' >> ../txt/C_capseq_TSS_capCounts_and_unique_snatch_counts_combined_time.txt
