for segment in {"1","2","3","4","5","6","7","8"}
do
	for time in {"A","B","C","D","E","F","G","H"}
	do
  	     	totalReads=`samtools view -F 4 ../bowtie_out/112213flucap_${time}.sorted.bam | wc -l`
		vRNA=`samtools view -f 16 ../bowtie_out/112213flucap_${time}.sorted.bam seg${segment}_B59FULL | wc -l`
		mRNA=`samtools view -F 16 ../bowtie_out/112213flucap_${time}.sorted.bam seg${segment}_B59FULL | wc -l`
		echo -e "${segment}\\t${time}\\t${mRNA}\\t${vRNA}\\t${totalReads}\\t" >> ../txt/reads_vRNA_mRNA.txt
	done
	echo "" >> ../txt/reads_vRNA_mRNA.txt
done
for time in {"6h","12h","24h","48h","mock12h"}
do
	echo ${time} >> ../txt/tpms_vRNA_mRNA_RNAseq.txt
	grep "^seg" ../stringtie_out/stringtie_out_${time}.gtf | awk -F"\t" '{ if ($3 == "transcript" && $7 == "+") {split($NF,arr," "); if (arr[length(arr)-1] == "TPM") {print $1"\t"arr[length(arr)]}}}' | sed 's/\"//g' | sed 's/;//g' | sort -k1,1 | bedtools groupby -g 1 -c 2 -o sum >> ../txt/tpms_vRNA_mRNA_RNAseq.txt
	grep "^seg" ../stringtie_out/stringtie_out_${time}.gtf | awk -F"\t" '{ if ($3 == "transcript" && $7 == "-") {split($NF,arr," "); if (arr[length(arr)-1] == "TPM") {print $1"\t"arr[length(arr)]}}}' | sed 's/\"//g' | sed 's/;//g' | sort -k1,1 | bedtools groupby -g 1 -c 2 -o sum >> ../txt/tpms_vRNA_mRNA_RNAseq.txt
done
