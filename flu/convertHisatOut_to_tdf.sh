for time in {"6h","12h","24h","48h","mock12h"}
do
echo """
module load IGVTools/2.3.31
igvtools count -z 5 -w 10 -e 200 /home/si14w/gnearline/flu/hisat2_out_${time}.sorted.bam   /home/si14w/gnearline/flu/hisat2_out/hisat2_out_${time}.sorted.tdf ../txt/hg19_plus_flu.chrom.sizes
""" > /home/si14w/gnearline/flu/bsubFiles/bamToTdf.${time}.bsub
done
