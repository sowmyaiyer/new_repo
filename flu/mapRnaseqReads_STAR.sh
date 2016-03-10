#cat biocore_ucsc.gtf flu.gtf > /home/si14w/gnearline/flu/txt/biocore_ucsc_plus_flu.gtf
# Dont really need to convert to bed, doing it just in case
#Rscript getGenesFromGtf.R /home/si14w/gnearline/flu/txt/biocore_ucsc_plus_flu.gtf /home/si14w/gnearline/flu/txt/biocore_ucsc_plus_flu.bed
#sort -k1,1 -k2,2n /home/si14w/gnearline/flu/txt/biocore_ucsc_plus_flu.bed > /home/si14w/gnearline/flu/txt/biocore_ucsc_plus_flu.sorted.bed
## Call STAR and cufflinks

#mkdir /project/umw_flustore_schiffer/STAR_out
for time in {"6h","12h","24h","48h","mock12h"}
do
echo """
module load star/2.4.2a
STAR --runThreadN 8 --genomeDir /project/umw_flustore_schiffer/STAR_index_hg19_and_fluWithSnp_combined --readFilesIn /farline/umw_flustore_schiffer/rnaseq/04SEP15.A549_B59_${time}.fastq.gz --readFilesCommand zcat --outFileNamePrefix /project/umw_flustore_schiffer/STAR_out/STAR_out_${time} --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outWigType bedGraph --outTmpDir /project/umw_flustore_schiffer/STAR_out/STAR_tmp_${time}
""" > ../bsubFiles/mapRnaseqReads_STAR_1st_pass.${time}.bsub
done

## 2nd pass mapping - multisample
## Collected sjdb from previous run, pass all of them to --sjdbFileChrStartEnd
for time in {"6h","12h","24h","48h","mock12h"}
do
echo """
module load star/2.4.2a
module load cufflinks/2.2.1
module load IGVTools/2.3.31

mkdir /home/si14w/gnearline/flu/cufflinks_out_${time}

STAR --runThreadN 8 --genomeDir /project/umw_flustore_schiffer/STAR_index_hg19_and_fluWithSnp_combined --readFilesIn /farline/umw_flustore_schiffer/rnaseq/04SEP15.A549_B59_${time}.fastq.gz --sjdbFileChrStartEnd /project/umw_flustore_schiffer/STAR_out/STAR_out_6hSJ.out.tab /project/umw_flustore_schiffer/STAR_out/STAR_out_12hSJ.out.tab /project/umw_flustore_schiffer/STAR_out/STAR_out_24hSJ.out.tab /project/umw_flustore_schiffer/STAR_out/STAR_out_48hSJ.out.tab /project/umw_flustore_schiffer/STAR_out/STAR_out_mock12hSJ.out.tab --readFilesCommand zcat --outFileNamePrefix /project/umw_flustore_schiffer/STAR_out/STAR_out_2ndpass_${time} --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outWigType bedGraph --outTmpDir /project/umw_flustore_schiffer/STAR_out/STAR_tmp_2ndpass_${time}
igvtools count -z 5 -w 25 -e 250 /project/umw_flustore_schiffer/STAR_out/STAR_out_2ndpass_${time}Aligned.sortedByCoord.out.bam  /project/umw_flustore_schiffer/STAR_out/STAR_out_2ndpass_${time}Aligned.sortedByCoord.out.tdf  hg19

cufflinks -p 8 -u -g $HOME/gnearline/flu/txt/biocore_ucsc_plus_flu.gtf --library-type fr-firststrand -o /home/si14w/gnearline/flu/cufflinks_out_${time} /project/umw_flustore_schiffer/STAR_out/STAR_out_2ndpass_${time}Aligned.sortedByCoord.out.bam
""" > ../bsubFiles/mapRnaseqReads_STAR_2ndpass_${time}.bsub
done
