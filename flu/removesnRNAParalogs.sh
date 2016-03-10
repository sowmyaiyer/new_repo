blastn -outfmt 6 -query ../txt/hg19_snRNA_and_snoRNA.sequences.fa -subject ../txt/hg19_snRNA_and_snoRNA.sequences.fa > ../txt/blast_output_snRNAs.txt
awk '{ if (($1 != $2) && ($3 == 100) && ($4 ==100) && ($5 == 0)) print }' ../txt/blast_output_snRNAs.txt | awk '{ print $2}' | sort | uniq  > ../txt/snRNAs_to_be_deleted.txt
grep "gl000" ../txt/hg19_snRNA_and_snoRNA.sequences.fa >> ../txt/snRNAs_to_be_deleted.txt
grep "coxx"  ../txt/hg19_snRNA_and_snoRNA.sequences.fa >> ../txt/snRNAs_to_be_deleted.txt
grep "mann"  ../txt/hg19_snRNA_and_snoRNA.sequences.fa >> ../txt/snRNAs_to_be_deleted.txt
grep "mcf"  ../txt/hg19_snRNA_and_snoRNA.sequences.fa >> ../txt/snRNAs_to_be_deleted.txt
grep "qbl"  ../txt/hg19_snRNA_and_snoRNA.sequences.fa >> ../txt/snRNAs_to_be_deleted.txt
grep "ssto"  ../txt/hg19_snRNA_and_snoRNA.sequences.fa >> ../txt/snRNAs_to_be_deleted.txt
sort -k1,1 ../txt/snRNAs_to_be_deleted.txt > ../txt/snRNAs_to_be_deleted.sorted.txt
awk '{ 
	if (index($0,">") == 1)  
		printf("%s\t",$0)
	else
		printf("%s\n",$0)
	}' ../txt/hg19_snRNA_and_snoRNA.sequences.fa | sed 's/>//g' | sort -k1,1 > ../txt/hg19_snRNA_and_snoRNA.sequences.columns.txt
join -v 1 ../txt/hg19_snRNA_and_snoRNA.sequences.columns.txt ../txt/snRNAs_to_be_deleted.sorted.txt -1 1 -2 1 | sed 's/ /\t/g' | awk '{ printf(">%s\n%s\n", $1,$2)}' > ../txt/hg19_snRNA_and_snoRNA.sequences.paralogs_removed.fa
