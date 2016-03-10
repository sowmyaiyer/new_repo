module load fastqc/0.10.1
for nm in `awk '{ print $1}' ../txt/chip01132016_sampleInfo.txt`
do
	filename=`grep -w $nm ../txt/chip01132016_sampleInfo.txt | awk '{ print $NF}'`
	echo $nm $filename
        cat /farline/umw_garberlab/human/DC/ChIP-Seq/ChIP_seq_1_13_2016/ChIP_test-27815805/$nm-3*/*R1_001.fastq.gz > /farline/umw_garberlab/human/DC/ChIP-Seq/ChIP_seq_1_13_2016/$filename.R1.fastq.gz
	cat /farline/umw_garberlab/human/DC/ChIP-Seq/ChIP_seq_1_13_2016/ChIP_test-27815805/$nm-3*/*R2_001.fastq.gz > /farline/umw_garberlab/human/DC/ChIP-Seq/ChIP_seq_1_13_2016/$filename.R2.fastq.gz
done
