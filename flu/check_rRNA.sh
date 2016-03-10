for time in {"A","B","C","D","E","F","G","H"}
do
        echo """
        module load bowtie2/2-2.1.0
        bowtie2 -x /share/data/umw_biocore/genome_data/human/hg19/commondb/rRNA/rRNA -U /project/umw_garberlab/narayanan/112213flucap/rawData_toFluAndtoHuman/${time}_Nocodes.fastq.gz -S ../bowtie_out/112213flucap_rRNA_${time}.sam
        samtools view -bS -o ../bowtie_out/112213flucap_rRNA_${time}.bam ../bowtie_out/112213flucap_rRNA_${time}.sam
        samtools sort -T ../bowtie_out/112213flucap_rRNA_${time}.sorted -o ../bowtie_out/112213flucap_rRNA_${time}.sorted.bam ../bowtie_out/112213flucap_rRNA_${time}.bam
        samtools index ../bowtie_out/112213flucap_rRNA_${time}.sorted.bam
        rm  ../bowtie_out/112213flucap_rRNA_${time}.sam  ../bowtie_out/112213flucap_rRNA_${time}.bam
        """ > ../bsubFiles/flu_capseq_bowtie_rRNA_${time}.bsub
done
