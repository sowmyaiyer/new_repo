for time in {"6h","12h","24h","48h","mock12h"}
do
echo """
module load IGVTools/2.3.31
#samtools index /project/umw_flustore_schiffer/STAR_out/STAR_out_${time}Aligned.sortedByCoord.out.bam
igvtools count -z 5 -w 25 -e 250 /project/umw_flustore_schiffer/STAR_out/STAR_out_${time}Aligned.sortedByCoord.out.bam  /project/umw_flustore_schiffer/STAR_out/STAR_out_${time}Aligned.sortedByCoord.out.tdf ../txt/hg19_plus_flu.chrom.sizes """ > /home/si14w/gnearline/flu/bsubFiles/bamToTdf.${time}.bsub
done
