for fastqgz in /farline/umw_garberlab/human/DC/ChIP-Seq/Chip_seq_12_15_15/Unindexed\ Reads-27272248/Undetermined\ from\ 151215_NB501205_0011_AHL5VYBGXX-31584787/*.fastq.gz
do
	echo $fastqgz
	fname=`basename "$fastqgz"`
echo """zcat \"${fastqgz}\" | awk '{ 
	if (index(\$0,\"@NB\") == 0) 
		printf(\"\\t%s\",\$0)
	else if (NR > 1)
		printf(\"\\n%s\",\$0)
	else if (NR == 1)
		printf(\"%s\",\$0)
	} 
	END{ printf(\"\\n\") } ' | awk -vnm=${fname} -F\"\\t\" '{ 
		split(\$1,readname,\":\"); 
		ind=substr(readname[length(readname)],1,6);
		print \$0 >> \"/home/si14w/gfarline/delete/\"nm\".\"ind\".txt\"
	}' 
""" > ../bsubFiles/demultiplex.$fname.bsub
done
#while read line
#do
#	bc=`echo ${line} | awk '{ print $1}'`
#	sample=`echo ${line} | awk '{ print $2}'`
#	echo ${bc}
#	sed 's/\t/\n/g' /home/si14w/gfarline/delete/${bc}.R1.txt | gzip > /farline/umw_garberlab/human/DC/ChIP-Seq/Chip_seq_12_15_15/demultiplexed/${sample}.R1.fastq.gz
#done < barcodes_121515.txt	
