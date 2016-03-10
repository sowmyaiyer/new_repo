for time in {"A","B","C","D","E","F","G","H"}
do
	echo $time
	sort -k4,4 /home/si14w/gnearline/flu/txt/112213flucap_${time}.readCounts_150bp.rpm.allRNA.genes.bed | awk '{ print $NF}' > /home/si14w/gnearline/flu/txt/112213flucap_${time}.readCountsOnly.genes.txt
done
sort -k4,4 /home/si14w/gnearline/flu/txt/112213flucap_${time}.readCounts_150bp.rpm.allRNA.genes.bed | awk '{ print $4}' > ../txt/geneList.genes.txt
echo gene 6h 6h_control 12h 12h_control 24h 24h_control 48h 48h_control > ../txt/capseq_rpm_allRNA.genes.txt
paste ../txt/geneList.genes.txt /home/si14w/gnearline/flu/txt/112213flucap_A.readCountsOnly.genes.txt /home/si14w/gnearline/flu/txt/112213flucap_E.readCountsOnly.genes.txt /home/si14w/gnearline/flu/txt/112213flucap_B.readCountsOnly.genes.txt /home/si14w/gnearline/flu/txt/112213flucap_F.readCountsOnly.genes.txt /home/si14w/gnearline/flu/txt/112213flucap_C.readCountsOnly.genes.txt /home/si14w/gnearline/flu/txt/112213flucap_G.readCountsOnly.genes.txt /home/si14w/gnearline/flu/txt/112213flucap_D.readCountsOnly.genes.txt /home/si14w/gnearline/flu/txt/112213flucap_H.readCountsOnly.genes.txt >>  ../txt/capseq_rpm_allRNA.genes.txt
