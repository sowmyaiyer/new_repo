for time in {"A","B","C","D","E","F","G","H"}
do
echo """Rscript fluSnatch.R /home/si14w/gnearline/flu/txt/${time}_capseq_TSS_capCounts_and_unique_snatch_counts.txt  /home/si14w/gnearline/flu/txt/112213flucap_${time}.multi_mapping_read_counts.proper.txt ${time} """ > ../bsubFiles/fluSnatch.capseqTSS.${time}.bsub
done
