awk '{ print $7}' ../data/H7hESC.DNase.Stam.narrowPeak > ../txt/DNase_scores_all.txt
Rscript getAtacScorePercentiles.R ../txt/DNase_scores_all.txt ../txt/dnase_quantile_values.txt
echo "bin numbers total" > dnase_scores_percentiles.txt
for n in {1..100}
do
	echo $n
	greater_than=`head -$((n+1)) ../txt//dnase_quantile_values.txt | tail -1`
	less_than=`head -$n ../txt/dnase_quantile_values.txt | tail -1`
	echo $greater_than $less_than 
	total=`awk \-vgreater_than=$greater_than \-vless_than=$less_than '{ if ($7 > greater_than && $7 < less_than) print }' ../data/H7hESC.DNase.Stam.narrowPeak | wc -l`
	number_intersecting=`awk \-vgreater_than=$greater_than \-vless_than=$less_than '{ if ($7 > greater_than && $7 < less_than) print }' ../data/H7hESC.DNase.Stam.narrowPeak | bedtools intersect -a stdin -b ../data/H1hesc.H3K27ac.ENCODE.bed -u | wc -l`
	echo $greater_than $number_intersecting $total >> dnase_scores_percentiles.txt
done


bedtools intersect -a ../data/H7hESC.DNase.Stam.narrowPeak -b ../data/GencodeV19.promoter.bed -u | bedtools intersect -a stdin -b ../data/H1hesc.H3K27ac.ENCODE.bed -u | awk '{ print $NF}' > ../txt/DNase_scores_intersecting_H3K27ac.proximal.txt
bedtools intersect -a ../data/H7hESC.DNase.Stam.narrowPeak -b ../data/GencodeV19.promoter.bed -v | bedtools intersect -a stdin -b ../data/H1hesc.H3K27ac.ENCODE.bed -u | awk '{ print $NF}' > ../txt/DNase_scores_intersecting_H3K27ac.distal.txt

bedtools intersect -a ../data/H7hESC.DNase.Stam.narrowPeak -b ../data/GencodeV19.promoter.bed -u | bedtools intersect -a stdin -b ../data/H1hesc.H3K27ac.ENCODE.bed -v | awk '{ print $NF}' > ../txt/DNase_scores_nonintersecting_H3K27ac.proximal.txt
bedtools intersect -a ../data/H7hESC.DNase.Stam.narrowPeak -b ../data/GencodeV19.promoter.bed -v | bedtools intersect -a stdin -b ../data/H1hesc.H3K27ac.ENCODE.bed -v | awk '{ print $NF}' > ../txt/DNase_scores_nonintersecting_H3K27ac.distal.txt

Rscript DNase_scores_proximal_distal.R
