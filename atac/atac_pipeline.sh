atac_peak_file=$1
#atac_peak_file=/project/umw_garberlab/tabakb/Maehr_AtacSeq/HumanES/22JUL14_PE50_C4MU1AC-B/Sample_HI_ATAC-Seq/Alignments/macs/HI_ATAC-Seq_NoIndex_L002.Merged.cutadapt.sorted_peaks.bed
output_dir_name=$2
echo $output_dir_name
mkdir ../results/${output_dir_name}
bedtools slop -i ../data/gencodeV19_TSS.bed -l 1000 -r 500 -s -g ~/TR/hg19.genome > ../data/GencodeV19.promoter.bed
bedtools intersect -a ${atac_peak_file} -b  ../data/GencodeV19.promoter.bed -v | sort -k1,1 -k2,2n > ../data/H1_ATAC.distal.sorted.bed

cat ../data/H1hesc.H3K27ac.ENCODE.bed ../data/H7_H3K27ac.bed | sort -k1,1 -k2,2n  | awk '{ print $1"\t"$2"\t"$3 }' | bedtools merge -i stdin > ../data/H1_and_H7_H3K27ac_merged.bed
bedtools intersect -a ../data/H1_ATAC.distal.sorted.bed -b ../data/H1_and_H7_H3K27ac_merged.bed -wao  | bedtools groupby -i stdin -c 9 -o max | awk 'BEGIN{print "label"}{ if  ($NF == 0) label=0; else label=1; print label}' > ../txt/enhancer_distal_positive_negative_labels.txt

awk '{ print $1"\t"$2"\t"$3"\t"$4}' ../data/H1_ATAC.distal.sorted.bed | bigWigAverageOverBed ~/farline/hot_mar2012/hg19.100way.phyloP100way.bw stdin ~/nearline/RECYCLE_BIN/out.tab -bedOut=../data/H1_ATAC.distal.cons.bed
awk '{ print $1"\t"$2"\t"$3"\t"$4}' ../data/H1_ATAC.distal.sorted.bed | bigWigAverageOverBed ../data/H1hesc.H3K9ac.bw stdin ~/nearline/RECYCLE_BIN/out.tab -bedOut=../data/H1_ATAC.distal.H3K9ac.bed
awk '{ print $1"\t"$2"\t"$3"\t"$4}' ../data/H1_ATAC.distal.sorted.bed | bigWigAverageOverBed ../data/H1hesc.H3K9me3.bw stdin ~/nearline/RECYCLE_BIN/out.tab -bedOut=../data/H1_ATAC.distal.H3K9me3.bed
awk '{ print $1"\t"$2"\t"$3"\t"$4}' ../data/H1_ATAC.distal.sorted.bed | bigWigAverageOverBed ../data/H1hesc.H3K4me1.bw stdin ~/nearline/RECYCLE_BIN/out.tab -bedOut=../data/H1_ATAC.distal.H3K4me1.bed
awk '{ print $1"\t"$2"\t"$3"\t"$4}' ../data/H1_ATAC.distal.sorted.bed | bigWigAverageOverBed ../data/H1hesc.H3K4me3.bw stdin ~/nearline/RECYCLE_BIN/out.tab -bedOut=../data/H1_ATAC.distal.H3K4me3.bed
awk '{ print $1"\t"$2"\t"$3"\t"$4}' ../data/H1_ATAC.distal.sorted.bed | bigWigAverageOverBed ../data/H1hesc.H3K27me3.bw stdin ~/nearline/RECYCLE_BIN/out.tab -bedOut=../data/H1_ATAC.distal.H3K27me3.bed
awk '{ print $1"\t"$2"\t"$3"\t"$4}' ../data/H1_ATAC.distal.sorted.bed | bigWigAverageOverBed ../data/H1hesc.H3K27ac.bw stdin ~/nearline/RECYCLE_BIN/out.tab -bedOut=../data/H1_ATAC.distal.H3K27ac.bed
awk '{ print $1"\t"$2"\t"$3"\t"$4}' ../data/H1_ATAC.distal.sorted.bed | bigWigAverageOverBed ../data/H1hesc.H3K4me2.bw stdin ~/nearline/RECYCLE_BIN/out.tab -bedOut=../data/H1_ATAC.distal.H3K4me2.bed

bedtools nuc -fi ~/my/hg19.fa -bed ../data/H1_ATAC.distal.sorted.bed | awk 'BEGIN{print ("gc")}{ if(NR > 1) print $7}' > ../txt/H1_ATAC.distal.GC.txt &

awk 'BEGIN{ print ("atac_score" ) }{ print $5 }' ../data/H1_ATAC.distal.sorted.bed >../txt/H1_ATAC.distal.atac_score.txt
sort -k1,1 -k2,2n ../data/H1_ATAC.distal.cons.bed | awk 'BEGIN{ print ("phylop_cons" ) }{ print $NF}'  >  ../txt/H1_ATAC.distal.phylop_cons.txt
bedtools closest -a ../data/H1_ATAC.distal.sorted.bed -b ../data/gencodeV19_TSS.bed -d -t "first" > ../txt/H1_ATAC_closest_tss.bed
awk 'BEGIN{ print ("dist_to_gene")} { print $NF}' ../txt/H1_ATAC_closest_tss.bed  > ../txt/H1_ATAC.distal.dist_to_gene.txt

sort -k1,1 -k2,2n ../data/H1_ATAC.distal.H3K4me1.bed | awk 'BEGIN{ print ("H3K4me1")} { print $NF}' > ../txt/H1_ATAC.distal.H3K4me1.txt
sort -k1,1 -k2,2n ../data/H1_ATAC.distal.H3K4me3.bed | awk 'BEGIN{ print ("H3K4me3")} { print $NF}' > ../txt/H1_ATAC.distal.H3K4me3.txt
sort -k1,1 -k2,2n ../data/H1_ATAC.distal.H3K9me3.bed | awk 'BEGIN{ print ("H3K9me3")} { print $NF}' > ../txt/H1_ATAC.distal.H3K9me3.txt
sort -k1,1 -k2,2n ../data/H1_ATAC.distal.H3K27me3.bed | awk 'BEGIN{ print ("H3K27me3")} { print $NF}' > ../txt/H1_ATAC.distal.H3K27me3.txt
sort -k1,1 -k2,2n ../data/H1_ATAC.distal.H3K9ac.bed | awk 'BEGIN{ print ("H3K9ac")} { print $NF}' > ../txt/H1_ATAC.distal.H3K9ac.txt
sort -k1,1 -k2,2n ../data/H1_ATAC.distal.H3K27ac.bed | awk 'BEGIN{ print ("H3K27ac")} { print $NF}' > ../txt/H1_ATAC.distal.H3K27ac.txt
sort -k1,1 -k2,2n ../data/H1_ATAC.distal.H3K4me2.bed | awk 'BEGIN{ print ("H3K4me2")} { print $NF}' > ../txt/H1_ATAC.distal.H3K4me2.txt

awk 'BEGIN{ print ("id")} { print $4 }' ../data/H1_ATAC.distal.sorted.bed > ../txt/H1_ATAC_seqids.txt


paste ../txt/H1_ATAC_seqids.txt ../txt/H1_ATAC.distal.atac_score.txt ../txt/H1_ATAC.distal.phylop_cons.txt ../txt/H1_ATAC.distal.dist_to_gene.txt  ../txt/H1_ATAC.distal.GC.txt ../txt/H1_ATAC.distal.H3K4me1.txt ../txt/H1_ATAC.distal.H3K4me3.txt ../txt/H1_ATAC.distal.H3K9me3.txt ../txt/H1_ATAC.distal.H3K27me3.txt ../txt/H1_ATAC.distal.H3K9ac.txt ../txt/H1_ATAC.distal.H3K27ac.txt ../txt/H1_ATAC.distal.H3K4me2.txt ../txt/enhancer_distal_positive_negative_labels.txt |  
awk 'BEGIN {
	OFS="\t"; 
	print "id\tatac_score\tphylop_cons\tdist_to_gene\tgc\tH3K4me1\tH3K4me3\tH3K9me3\tH3K27me3\tH3K9ac\tH3K27ac\tH3K4me2\tlabel"
	} 
	{ 
		if(NR > 1) { print }
	}' > ../results/${output_dir_name}/H1_ATAC_trainingset_withhismods.txt

Rscript atac_feature_comparison.R ${output_dir_name}
awk '{ print $1":"$2"-"$3}' ../data/H1_ATAC.distal.sorted.bed > ../txt/H1_ATAC_distal.seqList
twoBitToFa ~/my/hg19.2bit ../txt/H1_ATAC_peak_distal_sequences.fa -seqList=../txt/H1_ATAC_distal.seqList


awk '{
                if (substr($1,1,1) != ">")
                {
                        printf("%s",$0)
                } else if (NR != 1) {
                        printf("\n")
                }
            } END { printf("\n") }' ../txt/H1_ATAC_peak_distal_sequences.fa > ../results/${output_dir_name}/H1_ATAC_peak_distal_sequence.only.txt


#### TF analysis ######

if [[ -f ../results/${output_dir_name}/H1hesc_TF_intersection_stats.txt ]]; then
	rm ../results/${output_dir_name}/H1hesc_TF_intersection_stats.txt
fi
bedtools intersect -a ${atac_peak_file} -b ../data/H1_and_H7_H3K27ac_merged.bed -u > ../data/H1_ATAC_intersect_H3K27ac.bed
bedtools intersect -a ${atac_peak_file} -b ../data/H1_and_H7_H3K27ac_merged.bed -v > ../data/H1_ATAC_non_intersect_H3K27ac.bed
echo tfname peaks_intersecting_atac_enhancers peaks_intersecting_atac_nonenhancers peaks_intersecting_atac_both_enh_and_nonenh peaks_not_intersecting_atac total_peaks > ../results/${output_dir_name}/H1hesc_TF_intersection_stats.txt
rm ../data/H1hesc.AllTFPeaks.bed
for tfgz in ../data/ENCS*gz
do
	tfname=`basename $tfgz | cut -d"_" -f4 | cut -d"." -f1`
	echo $tfname
	zcat $tfgz | awk -vtfname=$tfname 'BEGIN{OFS="\t"}{ print $0,tfname}' >> ../data/H1hesc.AllTFPeaks.bed
	total_peaks=`zcat $tfgz | wc -l`
	zcat $tfgz | bedtools intersect -a stdin -b ../data/H1_ATAC_intersect_H3K27ac.bed -wa -u > ../data/$tfname.peaks_intersecting_atac_enhancers.bed
	zcat $tfgz | bedtools intersect -a stdin -b ../data/H1_ATAC_non_intersect_H3K27ac.bed -wa -u > ../data/$tfname.peaks_intersecting_atac_nonenhancers.bed
	peaks_intersecting_atac_both_enh_and_nonenh=`bedtools intersect -a ../data/$tfname.peaks_intersecting_atac_enhancers.bed -b ../data/$tfname.peaks_intersecting_atac_nonenhancers.bed -f 1.0 -r | wc -l`
	
	peaks_intersecting_atac_enhancers=`wc -l ../data/$tfname.peaks_intersecting_atac_enhancers.bed | cut -d" " -f1`
	peaks_intersecting_atac_nonenhancers=`wc -l ../data/$tfname.peaks_intersecting_atac_nonenhancers.bed | cut -d" " -f1`
	peaks_intersecting_atac_both_enh_and_nonenh=`bedtools intersect -a ../data/$tfname.peaks_intersecting_atac_enhancers.bed -b ../data/$tfname.peaks_intersecting_atac_nonenhancers.bed -f 1.0 -r | wc -l`
#	peaks_intersecting_atac_all=`zcat $tfgz | bedtools intersect -a stdin -b ../data/H1_ATAC.distal.sorted.bed -u | wc -l`
	peaks_not_intersecting_atac=`zcat $tfgz | bedtools intersect -a stdin -b ../data/H1_ATAC.distal.sorted.bed -v | wc -l`
	echo $tfname $peaks_intersecting_atac_enhancers $peaks_intersecting_atac_nonenhancers $peaks_intersecting_atac_both_enh_and_nonenh $peaks_not_intersecting_atac $total_peaks >> ../results/${output_dir_name}/H1hesc_TF_intersection_stats.txt
done
sed 's/ /\t/g' ../results/${output_dir_name}/H1hesc_TF_intersection_stats.txt > ../txt/tmp
mv ../txt/tmp ../results/${output_dir_name}/H1hesc_TF_intersection_stats.txt


# Now asking the reverse question. How many non-enhancer ATAC peaks are bound by each TF? How many enhancer ATAC peaks are bound by each TF?

total_enhancers=`wc -l ../data/H1_ATAC_intersect_H3K27ac.bed | cut -d" " -f1`
total_nonenhancers=`wc -l ../data/H1_ATAC_non_intersect_H3K27ac.bed | cut -d" " -f1`

echo  TF enhancers_bound_to_tf total_enhancers totalTFPeaks > ../results/${output_dir_name}/enhancer_intersecting_tf.txt
echo  TF nonenhancers_bound_to_tf total_nonenhancers totalTFPeaks > ../results/${output_dir_name}/nonenhancer_intersecting_tf.txt
for tfgz in ../data/ENCS*gz
do
	tfname=`basename $tfgz | cut -d"_" -f4 | cut -d"." -f1`
	totalPeaks=`zcat $tfgz| wc -l | cut -d" " -f1`
	enhancers_intersecting_tf=`zcat $tfgz | bedtools intersect -a ../data/H1_ATAC_intersect_H3K27ac.bed -b stdin -u | wc -l`
	nonenhancers_intersecting_tf=`zcat $tfgz | bedtools intersect -a ../data//H1_ATAC_non_intersect_H3K27ac.bed -b stdin -u | wc -l`
	echo $tfname $enhancers_intersecting_tf $total_enhancers $totalPeaks >> ../results/${output_dir_name}/enhancer_intersecting_tf.txt
	echo $tfname $nonenhancers_intersecting_tf $total_nonenhancers $totalPeaks >> ../results/${output_dir_name}/nonenhancer_intersecting_tf.txt
done

bedtools intersect -a ../data/H1_ATAC.distal.sorted.bed -b ../data/H1hesc.AllTFPeaks.bed -wa -wb |  awk 'BEGIN{OFS="\t"}{ print $16,$5}' > ../results/${output_dir_name}/TFs_and_ATAC_scores.txt


bedtools intersect -a ../data/H1hesc.AllTFPeaks.bed -b ../data/H1_ATAC_intersect_H3K27ac.bed | awk 'BEGIN{OFS="\t"}{ print $11,$7,"enhancer_ATAC"}' > ../results/${output_dir_name}/H1_TFPeaks_in_enhancer_atac_peaks.txt
bedtools intersect -a ../data/H1hesc.AllTFPeaks.bed -b ../data/H1_ATAC_non_intersect_H3K27ac.bed | awk 'BEGIN{OFS="\t"}{ print $11,$7,"non_enhancer_ATAC" }' > ../results/${output_dir_name}/H1_TFPeaks_in_nonenhancer_atac_peaks.txt
bedtools intersect -a ../data/H1hesc.AllTFPeaks.bed -b ../data/H1_ATAC_intersect_H3K27ac.bed -v | awk 'BEGIN{OFS="\t"}{ print $11,$7,"non_ATAC"}' > ../results/${output_dir_name}/H1_TFPeaks_not_in_atac_peaks.txt
Rscript tfPeakAnalysis.R ${output_dir_name}

## END TF analysis ###

## BEGIN expression analysis ###

bedtools intersect -a ${atac_peak_file}  -b ../data/H1_and_H7_H3K27ac_merged.bed -u > ../data/H1_ATAC_intersect_H3K27ac.bed
bedtools intersect -a ../data/H1_ATAC.distal.sorted.bed -b ../data/H1_and_H7_H3K27ac_merged.bed -v > ../data/H1_ATAC_non_intersect_H3K27ac.bed
cat ../data/GencodeV19.promoter.bed ../data/H1_ATAC.distal.sorted.bed | awk 'BEGIN{OFS="\t"}{ print $1,$2,$3}' > ../txt/excludable_ATAC_and_promoters.bed

awk '{ if ($3 == "transcript") print }' ../data/H1hesc.transcript.quant.gtf  | awk -f extractFpkmWithTss.awk > ../data/H1hesc.transcript_only.tss_and_fpkm.bed

bedtools closest -a  ../data/H1_ATAC_intersect_H3K27ac.bed -b ../data/H1hesc.transcript_only.tss_and_fpkm.bed | awk 'BEGIN{OFS="\t"}{print "enhancer",$9}' > ../results/${output_dir_name}/fpkm_genes_closest_to_atac_enhancers.txt
bedtools closest -a  ../data/H1_ATAC_non_intersect_H3K27ac.bed -b ../data/H1hesc.transcript_only.tss_and_fpkm.bed | awk 'BEGIN{OFS="\t"}{print "nonenhancer",$9}' > ../results/${output_dir_name}/fpkm_genes_closest_to_atac_nonenhancers.txt
bedtools shuffle -i ../data/H1_ATAC.distal.sorted.bed -g ../data/hg19.genome -excl ../txt/excludable_ATAC_and_promoters.bed | bedtools closest -a stdin -b ../data/H1hesc.transcript_only.tss_and_fpkm.bed -d | awk 'BEGIN{OFS="\t"}{if ($7 != -1) { print "random",$9} }' > ../results/${output_dir_name}/fpkm_genes_closest_to_random_distal_regions.txt


bedtools closest -a  ../data/H1_ATAC_intersect_H3K27ac.bed -b ../data/H1hesc.transcript_only.tss_and_fpkm.bed -d | awk 'BEGIN{OFS="\t"}{if ($10 <= 10000) { print "enhancer",$9} }' > ../results/${output_dir_name}/fpkm_genes_closest_to_atac_enhancers_within10k.txt
bedtools closest -a  ../data/H1_ATAC_non_intersect_H3K27ac.bed -b ../data/H1hesc.transcript_only.tss_and_fpkm.bed -d | awk 'BEGIN{OFS="\t"}{if ($10 <= 10000) { print "nonenhancer",$9} }' > ../results/${output_dir_name}/fpkm_genes_closest_to_atac_nonenhancers_within10k.txt
bedtools shuffle -i ../data/H1_ATAC.distal.sorted.bed -g ../data/hg19.genome -excl ../txt/excludable_ATAC_and_promoters.bed | bedtools closest -a stdin -b ../data/H1hesc.transcript_only.tss_and_fpkm.bed -d | awk 'BEGIN{OFS="\t"}{if ($10 <= 10000 && $7 != -1) { print "random",$9} }' > ../results/${output_dir_name}/fpkm_genes_closest_to_random_distal_regions_within10k.txt


bedtools closest -a  ../data/H1_ATAC_intersect_H3K27ac.bed -b ../data/H1hesc.transcript_only.tss_and_fpkm.bed -d | awk 'BEGIN{OFS="\t"}{if ($10 <= 5000) { print "enhancer",$9} }' > ../results/${output_dir_name}/fpkm_genes_closest_to_atac_enhancers_within5k.txt
bedtools closest -a  ../data/H1_ATAC_non_intersect_H3K27ac.bed -b ../data/H1hesc.transcript_only.tss_and_fpkm.bed -d | awk 'BEGIN{OFS="\t"}{if ($10 <= 5000) { print "nonenhancer",$9} }' > ../results/${output_dir_name}/fpkm_genes_closest_to_atac_nonenhancers_within5k.txt
bedtools shuffle -i ../data/H1_ATAC.distal.sorted.bed -g ../data/hg19.genome -excl ../txt/excludable_ATAC_and_promoters.bed | bedtools closest -a stdin -b ../data/H1hesc.transcript_only.tss_and_fpkm.bed -d | awk 'BEGIN{OFS="\t"}{if ($10 <= 5000 && $7 != -1) { print "random",$9} }' > ../results/${output_dir_name}/fpkm_genes_closest_to_random_distal_regions_within5k.txt
Rscript expressionAnalysis.R ${output_dir_name}
## END expression analysis ###
