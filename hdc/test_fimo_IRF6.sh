for chr in `seq 1 22`
do
	echo """ 
	fimo --text --thresh 0.00001 ../txt/M1266_1.02_IRF6.meme /share/data/umw_biocore/genome_data/human/hg19/chroms/chr${chr}.fa | awk 'BEGIN{OFS=\"\t\"} \$7~/^[0-9.]+\$/ {print \$7 >> \"/farline/umw_garberlab/iyers/fimo_test/fimo_test_M1266_1.02_IRF6_chr${chr}.pvalues_only.txt\";print \$2,\$3,\$4,\$1,\$6,\$5 >> \"/farline/umw_garberlab/iyers/fimo_test/fimo_test_M1266_1.02_IRF6_chr${chr}.locations.txt\" }'
	""" > ../bsubFiles/fimo_test_M1266_1.02_IRF6_chr${chr}.bsub
done
echo """ fimo --text --thresh 0.00001 ../txt/M1266_1.02_IRF6.meme /share/data/umw_biocore/genome_data/human/hg19/chroms/chrX.fa | awk 'BEGIN{OFS=\"\t\"} \$7~/^[0-9.]+\$/ {print \$7 >> \"/farline/umw_garberlab/iyers/fimo_test/fimo_test_M1266_1.02_IRF6_chrX.pvalues_only.txt\";print \$2,\$3,\$4,\$5 >> \"/farline/umw_garberlab/iyers/fimo_test/fimo_test_M1266_1.02_IRF6_chrX.locations.txt\" }'
""" > ../bsubFiles/fimo_test_M1266_1.02_IRF6_chrX.bsub


for chrfile in `\ls /home/si14w/gfarline/fimo_test/fimo_test_M1266_1.02_IRF6_chr*.pvalues_only.txt | sort -k1,1`
do
	echo $chrfile
	cat $chrfile >> /farline/umw_garberlab/iyers/fimo_test/fimo_test_M1266_1.02_IRF6.p_values.txt
done

for chrfile in `\ls /home/si14w/gfarline/fimo_test/fimo_test_M1266_1.02_IRF6_chr*.locations.txt | sort -k1,1`
do
        echo $chrfile
        cat $chrfile >> /farline/umw_garberlab/iyers/fimo_test/fimo_test_M1266_1.02_IRF6.locations.txt
done

qvalue --verbosity 4 /farline/umw_garberlab/iyers/fimo_test/fimo_test_M1266_1.02_IRF6.p_values.txt > /farline/umw_garberlab/iyers/fimo_test/fimo_test_M1266_1.02_IRF6_qvalue.txt
