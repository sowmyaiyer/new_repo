#Rscript generateRandomDNA.R > ../txt/random9_15mers.fa
#module load bowtie2/2-2.1.0
#bowtie2 -p 8 -x ../txt/TSS_plus_15.combined_time.fa -N 0 --norc -U ../txt/random9_15mers.fa -f --all -S ../txt/random10_15mer_aligned.sam
#grep -v "^@" ../txt/random10_15mer_aligned.sam | awk '{ print $1}' | sort | uniq -c | awk '{ print $1}' > ../txt/random_alignment_stats.txt
#Rscript plotRandomAlignmentStats.R ../txt/random_alignment_stats.txt ../results/random10_15mer.pdf

for time in {"A","B","C","D","E","F","G","H"}
do
	echo $time
	grep -v "^@" ../txt/${time}_capseq_TSS.combined_time.sam | awk '{ print $1}' | sort | uniq -c | awk '{ print $1}' > ../txt/${time}_alignment_stats.txt
	Rscript plotRandomAlignmentStats.R ../txt/${time}_alignment_stats.txt ../results/${time}_matches.pdf ${time}
done


