#sort -k1,1 -k2,2n ../data/H1_ATAC.distal.bed > ../data/H1_ATAC.distal.sorted.bed
zcat ../data/ENCSR000AMF_H1-hESC_ChIP-seq_CTCF.gz ../data/ENCSR000BLD_H1-hESC_ChIP-seq_RAD21.gz | sort -k1,1 -k2,2n | bedtools merge -i stdin | bedtools intersect -a ../data/H1_ATAC.distal.sorted.bed -b stdin -v | sort -k1,1 -k2,2n > ../data/H1_ATAC.distal_noCTCF.sorted.bed
awk '{ print $1"\t"$2"\t"$3"\t"$4}' ../data/H1_ATAC.distal_noCTCF.sorted.bed | bigWigAverageOverBed ~/farline/hot_mar2012/hg19.100way.phyloP100way.bw stdin ~/nearline/RECYCLE_BIN/out.tab -bedOut=/home/si14w/gnearline/ATAC_analysis/data/H1_ATAC.distal_noCTCF.cons.bed
awk '{ print $1"\t"$2"\t"$3"\t"$4}' ../data/H1_ATAC.distal_noCTCF.sorted.bed | bigWigAverageOverBed /home/si14w/gnearline/ATAC_analysis/data/H1hesc.H3K9ac.bw stdin ~/nearline/RECYCLE_BIN/out.tab -bedOut=/home/si14w/gnearline/ATAC_analysis/data/H1_ATAC.distal_noCTCF.H3K9ac.bed
awk '{ print $1"\t"$2"\t"$3"\t"$4}' ../data/H1_ATAC.distal_noCTCF.sorted.bed | bigWigAverageOverBed /home/si14w/gnearline/ATAC_analysis/data/H1hesc.H3K9me3.bw stdin ~/nearline/RECYCLE_BIN/out.tab -bedOut=/home/si14w/gnearline/ATAC_analysis/data/H1_ATAC.distal_noCTCF.H3K9me3.bed
awk '{ print $1"\t"$2"\t"$3"\t"$4}' ../data/H1_ATAC.distal_noCTCF.sorted.bed | bigWigAverageOverBed /home/si14w/gnearline/ATAC_analysis/data/H1hesc.H3K4me1.bw stdin ~/nearline/RECYCLE_BIN/out.tab -bedOut=/home/si14w/gnearline/ATAC_analysis/data/H1_ATAC.distal_noCTCF.H3K4me1.bed
awk '{ print $1"\t"$2"\t"$3"\t"$4}' ../data/H1_ATAC.distal_noCTCF.sorted.bed | bigWigAverageOverBed /home/si14w/gnearline/ATAC_analysis/data/H1hesc.H3K4me3.bw stdin ~/nearline/RECYCLE_BIN/out.tab -bedOut=/home/si14w/gnearline/ATAC_analysis/data/H1_ATAC.distal_noCTCF.H3K4me3.bed
awk '{ print $1"\t"$2"\t"$3"\t"$4}' ../data/H1_ATAC.distal_noCTCF.sorted.bed | bigWigAverageOverBed /home/si14w/gnearline/ATAC_analysis/data/H1hesc.H3K27me3.bw stdin ~/nearline/RECYCLE_BIN/out.tab -bedOut=/home/si14w/gnearline/ATAC_analysis/data/H1_ATAC.distal_noCTCF.H3K27me3.bed
awk '{ print $1"\t"$2"\t"$3"\t"$4}' ../data/H1_ATAC.distal_noCTCF.sorted.bed | bigWigAverageOverBed /home/si14w/gnearline/ATAC_analysis/data/H1hesc.H3K27ac.bw stdin ~/nearline/RECYCLE_BIN/out.tab -bedOut=/home/si14w/gnearline/ATAC_analysis/data/H1_ATAC.distal_noCTCF.H3K27ac.bed
awk '{ print $1"\t"$2"\t"$3"\t"$4}' ../data/H1_ATAC.distal_noCTCF.sorted.bed | bigWigAverageOverBed /home/si14w/gnearline/ATAC_analysis/data/H1hesc.H3K4me2.bw stdin ~/nearline/RECYCLE_BIN/out.tab -bedOut=/home/si14w/gnearline/ATAC_analysis/data/H1_ATAC.distal_noCTCF.H3K4me2.bed


#bedtools intersect -a /home/si14w/gnearline/ATAC_analysis/data/H1_ATAC.distal_noCTCF.sorted.bed -b ../data/H7_H3K27ac_distal.bed -wao  | bedtools groupby -i stdin -c 12 -o max | awk 'BEGIN{print "label"}{ if  ($NF == 0) label=0; else label=1; print label}' > ../txt/enhancer_distal_positive_negative_labels.txt
cat ../data/H1hesc.H3K27ac.ENCODE.bed ../data/H7_H3K27ac.bed | sort -k1,1 -k2,2n  | awk '{ print $1"\t"$2"\t"$3 }' | bedtools merge -i stdin > ../data/H1_and_H7_H3K27ac_merged.bed

bedtools intersect -a /home/si14w/gnearline/ATAC_analysis/data/H1_ATAC.distal_noCTCF.sorted.bed -b ../data/H1_and_H7_H3K27ac_merged.bed -wao  | bedtools groupby -i stdin -c 9 -o max | awk 'BEGIN{print "label"}{ if  ($NF == 0) label=0; else label=1; print label}' > ../txt/enhancer_distal_positive_negative_labels.txt
bedtools nuc -fi ~/my/hg19.fa -bed ../data/H1_ATAC.distal_noCTCF.sorted.bed | awk 'BEGIN{print ("seqLen")}{ if(NR > 1) print $14}' >  ../txt/H1_ATAC.distal_noCTCF.seqLen.txt &
bedtools nuc -fi ~/my/hg19.fa -bed ../data/H1_ATAC.distal_noCTCF.sorted.bed | awk 'BEGIN{print ("gc")}{ if(NR > 1) print $7}' > ../txt/H1_ATAC.distal_noCTCF.GC.txt &
bedtools nuc -fi ~/my/hg19.fa -bed ../data/H1_ATAC.distal_noCTCF.sorted.bed -pattern "CA" -C | awk 'BEGIN{ print ("ca" ) } { if (NR > 1) { printf("%.2f\n", $NF)} }' > ../txt/H1_ATAC.distal_noCTCF.CA.txt &
bedtools nuc -fi ~/my/hg19.fa -bed ../data/H1_ATAC.distal_noCTCF.sorted.bed -pattern "TG" -C | awk 'BEGIN{ print ("ca_rev" ) }{ if (NR > 1) { printf("%.2f\n", $NF)} }' > ../txt/H1_ATAC.distal_noCTCF.CA_rev.txt &
bedtools nuc -fi ~/my/hg19.fa -bed ../data/H1_ATAC.distal_noCTCF.sorted.bed -pattern "GA" -C | awk 'BEGIN{ print ("ga" ) } { if (NR > 1) { printf("%.2f\n", $NF)} }' > ../txt/H1_ATAC.distal_noCTCF.GA.txt &
bedtools nuc -fi ~/my/hg19.fa -bed ../data/H1_ATAC.distal_noCTCF.sorted.bed -pattern "TC" -C | awk 'BEGIN{ print ("ga_rev" ) } { if (NR > 1) { printf("%.2f\n", $NF)} }' > ../txt/H1_ATAC.distal_noCTCF.GA_rev.txt &
bedtools nuc -fi ~/my/hg19.fa -bed ../data/H1_ATAC.distal_noCTCF.sorted.bed -pattern "GC" -C | awk 'BEGIN{ print ("gc_dinuc" ) } { if (NR > 1) { printf("%.2f\n", $NF)} }' > ../txt/H1_ATAC.distal_noCTCF.GC_dinuc.txt &

awk 'BEGIN{ print ("atac_score" ) }{ print $5 }' ../data/H1_ATAC.distal_noCTCF.sorted.bed >../txt/H1_ATAC.distal_noCTCF.atac_score.txt
#awk 'BEGIN{ print ("atac_score" ) }{ print $5 }' ../txt/H7Hesc_ATAC_signal_in_ATAC_peaks.bed >../txt/H1_ATAC.distal_noCTCF.atac_score.txt
sort -k1,1 -k2,2n /home/si14w/gnearline/ATAC_analysis/data/H1_ATAC.distal_noCTCF.cons.bed | awk 'BEGIN{ print ("phylop_cons" ) }{ print $NF}'  >  ../txt/H1_ATAC.distal_noCTCF.phylop_cons.txt
bedtools closest -a ../data/H1_ATAC.distal_noCTCF.sorted.bed -b ../data/gencodeV19_TSS.bed -d -t "first" > ../txt/H1_ATAC_closest_tss.bed
awk 'BEGIN{ print ("dist_to_gene")} { print $NF}' ../txt/H1_ATAC_closest_tss.bed  > ../txt/H1_ATAC.distal_noCTCF.dist_to_gene.txt

sort -k1,1 -k2,2n /home/si14w/gnearline/ATAC_analysis/data/H1_ATAC.distal_noCTCF.H3K4me1.bed | awk 'BEGIN{ print ("H3K4me1")} { print $NF}' > ../txt/H1_ATAC.distal_noCTCF.H3K4me1.txt
sort -k1,1 -k2,2n /home/si14w/gnearline/ATAC_analysis/data/H1_ATAC.distal_noCTCF.H3K4me3.bed | awk 'BEGIN{ print ("H3K4me3")} { print $NF}' > ../txt/H1_ATAC.distal_noCTCF.H3K4me3.txt
sort -k1,1 -k2,2n /home/si14w/gnearline/ATAC_analysis/data/H1_ATAC.distal_noCTCF.H3K9me3.bed | awk 'BEGIN{ print ("H3K9me3")} { print $NF}' > ../txt/H1_ATAC.distal_noCTCF.H3K9me3.txt
sort -k1,1 -k2,2n /home/si14w/gnearline/ATAC_analysis/data/H1_ATAC.distal_noCTCF.H3K27me3.bed | awk 'BEGIN{ print ("H3K27me3")} { print $NF}' > ../txt/H1_ATAC.distal_noCTCF.H3K27me3.txt
sort -k1,1 -k2,2n /home/si14w/gnearline/ATAC_analysis/data/H1_ATAC.distal_noCTCF.H3K9ac.bed | awk 'BEGIN{ print ("H3K9ac")} { print $NF}' > ../txt/H1_ATAC.distal_noCTCF.H3K9ac.txt
sort -k1,1 -k2,2n /home/si14w/gnearline/ATAC_analysis/data/H1_ATAC.distal_noCTCF.H3K27ac.bed | awk 'BEGIN{ print ("H3K27ac")} { print $NF}' > ../txt/H1_ATAC.distal_noCTCF.H3K27ac.txt
sort -k1,1 -k2,2n /home/si14w/gnearline/ATAC_analysis/data/H1_ATAC.distal_noCTCF.H3K4me2.bed | awk 'BEGIN{ print ("H3K4me2")} { print $NF}' > ../txt/H1_ATAC.distal_noCTCF.H3K4me2.txt

awk 'BEGIN{ print ("id")} { print $4 }' ../data/H1_ATAC.distal_noCTCF.sorted.bed > ../txt/H1_ATAC_seqids.txt


paste ../txt/H1_ATAC_seqids.txt ../txt/H1_ATAC.distal_noCTCF.atac_score.txt ../txt/H1_ATAC.distal_noCTCF.phylop_cons.txt ../txt/H1_ATAC.distal_noCTCF.dist_to_gene.txt  ../txt/H1_ATAC.distal_noCTCF.GC.txt ../txt/H1_ATAC.distal_noCTCF.CA.txt ../txt/H1_ATAC.distal_noCTCF.CA_rev.txt ../txt/H1_ATAC.distal_noCTCF.GA.txt ../txt/H1_ATAC.distal_noCTCF.GA_rev.txt ../txt/H1_ATAC.distal_noCTCF.GC_dinuc.txt ../txt/H1_ATAC.distal_noCTCF.seqLen.txt ../txt/H1_ATAC.distal_noCTCF.H3K4me1.txt ../txt/H1_ATAC.distal_noCTCF.H3K4me3.txt ../txt/H1_ATAC.distal_noCTCF.H3K9me3.txt ../txt/H1_ATAC.distal_noCTCF.H3K27me3.txt ../txt/H1_ATAC.distal_noCTCF.H3K9ac.txt ../txt/H1_ATAC.distal_noCTCF.H3K27ac.txt ../txt/H1_ATAC.distal_noCTCF.H3K4me2.txt ../txt/enhancer_distal_positive_negative_labels.txt |  awk 'BEGIN{print "id\tatac_score\tphylop_cons\tdist_to_gene\tgc\tca\tga\tgc_dinuc\tH3K4me1\tH3K4me3\tH3K9me3\tH3K27me3\tH3K9ac\tH3K27ac\tH3K4me2\tlabel"} { if(NR > 1) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"($6+$7)/$11"\t"($8+$9)/$11"\t"$10/$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$NF} }' > ../txt/H1_ATAC.distal_noCTCF.trainingset_withhismods.txt

awk '{ print $1":"$2"-"$3}' ../data/H1_ATAC.distal_noCTCF.sorted.bed > ../txt/H1_ATAC_distal_noCTCF.seqList
twoBitToFa ~/my/hg19.2bit ../txt/H1_ATAC_peak_distal_noCTCF_sequences.fa -seqList=../txt/H1_ATAC_distal_noCTCF.seqList


awk '{
                if (substr($1,1,1) != ">")
                {
                        printf("%s",$0)
                } else if (NR != 1) {
                        printf("\n")
                }
            } END { printf("\n") }' ../txt/H1_ATAC_peak_distal_noCTCF_sequences.fa > ../txt/H1_ATAC_peak_distal_noCTCF_sequence.only.txt


