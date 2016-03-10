while read line
do
       bc=`echo ${line} | awk '{ print $1}'`
       sample=`echo ${line} | awk '{ print $2}'`
       echo ${bc}
       cat $HOME/gfarline/delete/Undetermined_S0_L*_R1_*.fastq.gz.${bc}.txt | sed 's/\t/\n/g' | gzip > /farline/umw_garberlab/human/DC/ChIP-Seq/Chip_seq_12_15_15/demultiplexed/${sample}.R1.fastq.gz
       cat $HOME/gfarline/delete/Undetermined_S0_L*_R2_*.fastq.gz.${bc}.txt | sed 's/\t/\n/g' | gzip > /farline/umw_garberlab/human/DC/ChIP-Seq/Chip_seq_12_15_15/demultiplexed/${sample}.R2.fastq.gz
done < barcodes_121515.txt
