#bedtools slop -i /project/umw_garberlab/tabakb/Maehr_AtacSeq/HumanES/22JUL14_PE50_C4MU1AC-B/Sample_HI_ATAC-Seq/Alignments/macs/HI_ATAC-Seq_NoIndex_L002.Merged.cutadapt.sorted_peaks.bed -b 50 -g ~/TR/hg19.genome >  H7_DNase.bed
#bedtools slop -i /home/si14w/nearline/enhancer_predictions/gencodeV19_TSS.bed -l 1000 -r 500 -s -g ~/TR/hg19.genome > GencodeV19.promoter.bed
bedtools intersect -a ../data/H7hESC.DNase.Stam.narrowPeak -b  ../data/GencodeV19.promoter.bed -v > ../data/H7_DNase.distal.bed
#bedtools intersect -a /project/umw_garberlab/tabakb/Maehr_AtacSeq/HumanES/22JUL14_PE50_C4MU1AC-B/Sample_HI_ATAC-Seq/Alignments/macs/HI_ATAC-Seq_NoIndex_L002.Merged.cutadapt.sorted_peaks.bed -b  GencodeV19.promoter.bed -v > H7_DNase.distal.bed

awk '{ print $1"\t"$2"\t"$3"\t"NR}' ../data/H7_DNase.distal.bed | bigWigAverageOverBed ~/farline/hot_mar2012/hg19.100way.phyloP100way.bw stdin ~/nearline/RECYCLE_BIN/out.tab -bedOut=/home/si14w/gnearline/ATAC_analysis/data/H7_DNase.distal.cons.bed
awk '{ print $1"\t"$2"\t"$3"\t"NR}' ../data/H7_DNase.distal.bed | bigWigAverageOverBed /home/si14w/gnearline/ATAC_analysis/data/H1hesc.H3K9ac.bw stdin ~/nearline/RECYCLE_BIN/out.tab -bedOut=/home/si14w/gnearline/ATAC_analysis/data/H7_DNase.distal.H3K9ac.bed
awk '{ print $1"\t"$2"\t"$3"\t"NR}' ../data/H7_DNase.distal.bed | bigWigAverageOverBed /home/si14w/gnearline/ATAC_analysis/data/H1hesc.H3K9me3.bw stdin ~/nearline/RECYCLE_BIN/out.tab -bedOut=/home/si14w/gnearline/ATAC_analysis/data/H7_DNase.distal.H3K9me3.bed
awk '{ print $1"\t"$2"\t"$3"\t"NR}' ../data/H7_DNase.distal.bed | bigWigAverageOverBed /home/si14w/gnearline/ATAC_analysis/data/H1hesc.H3K4me1.bw stdin ~/nearline/RECYCLE_BIN/out.tab -bedOut=/home/si14w/gnearline/ATAC_analysis/data/H7_DNase.distal.H3K4me1.bed
awk '{ print $1"\t"$2"\t"$3"\t"NR}' ../data/H7_DNase.distal.bed | bigWigAverageOverBed /home/si14w/gnearline/ATAC_analysis/data/H1hesc.H3K4me3.bw stdin ~/nearline/RECYCLE_BIN/out.tab -bedOut=/home/si14w/gnearline/ATAC_analysis/data/H7_DNase.distal.H3K4me3.bed
awk '{ print $1"\t"$2"\t"$3"\t"NR}' ../data/H7_DNase.distal.bed | bigWigAverageOverBed /home/si14w/gnearline/ATAC_analysis/data/H1hesc.H3K27me3.bw stdin ~/nearline/RECYCLE_BIN/out.tab -bedOut=/home/si14w/gnearline/ATAC_analysis/data/H7_DNase.distal.H3K27me3.bed

#bedtools intersect -a /home/si14w/gnearline/ATAC_analysis/data/H7_DNase.distal.bed -b ../data/H7_H3K27ac_distal.bed -wao  | bedtools groupby -i stdin -c 12 -o max | awk 'BEGIN{print "label"}{ if  ($NF == 0) label=0; else label=1; print label}' > ../txt/enhancer_distal_positive_negative_labels.txt
cat ../data/H1hesc.H3K27ac.ENCODE.distal.bed ../data/H7_H3K27ac_distal.bed | sort -k1,1 -k2,2n  | awk '{ print $1"\t"$2"\t"$3 }' | bedtools merge -i stdin > ../data/H1_and_H7_H3K27ac_merged.distal.bed

bedtools intersect -a /home/si14w/gnearline/ATAC_analysis/data/H7_DNase.distal.bed -b ../data/H1_and_H7_H3K27ac_merged.distal.bed -wao  | bedtools groupby -i stdin -c 12 -o max | awk 'BEGIN{print "label"}{ if  ($NF == 0) label=0; else label=1; print label}' > ../txt/enhancer_distal_positive_negative_labels.txt

bedtools nuc -fi ~/my/hg19.fa -bed ../data/H7_DNase.distal.bed | awk 'BEGIN{print ("seqLen")}{ if(NR > 1) print $14}' >  ../txt/H7_DNase.distal.seqLen.txt &
bedtools nuc -fi ~/my/hg19.fa -bed ../data/H7_DNase.distal.bed | awk 'BEGIN{print ("gc")}{ if(NR > 1) print $7}' > ../txt/H7_DNase.distal.GC.txt &

awk 'BEGIN{ print ("atac_score" ) }{ print $7 }' ../data/H7_DNase.distal.bed >../txt/H7_DNase.distal.atac_score.txt
#awk 'BEGIN{ print ("atac_score" ) }{ print $5 }' ../txt/H7Hesc_ATAC_signal_in_ATAC_peaks.bed >../txt/H7_DNase.distal.atac_score.txt
awk 'BEGIN{ print ("phylop_cons" ) }{ print $NF}' /home/si14w/gnearline/ATAC_analysis/data/H7_DNase.distal.cons.bed >  ../txt/H7_DNase.distal.phylop_cons.txt
bedtools closest -a ../data/H7_DNase.distal.bed -b /home/si14w/nearline/enhancer_predictions/gencodeV19_TSS.bed -d -t "first" > ../txt/H7_DNase_closest_tss.bed
awk 'BEGIN{ print ("dist_to_gene")} { print $NF}' ../txt/H7_DNase_closest_tss.bed  > ../txt/H7_DNase.distal.dist_to_gene.txt

awk 'BEGIN{ print ("H3K4me1")} { print $NF}' /home/si14w/gnearline/ATAC_analysis/data/H7_DNase.distal.H3K4me1.bed  > ../txt/H7_DNase.distal.H3K4me1.txt
awk 'BEGIN{ print ("H3K4me3")} { print $NF}' /home/si14w/gnearline/ATAC_analysis/data/H7_DNase.distal.H3K4me3.bed  > ../txt/H7_DNase.distal.H3K4me3.txt
awk 'BEGIN{ print ("H3K9me3")} { print $NF}' /home/si14w/gnearline/ATAC_analysis/data/H7_DNase.distal.H3K9me3.bed  > ../txt/H7_DNase.distal.H3K9me3.txt
awk 'BEGIN{ print ("H3K27me3")} { print $NF}' /home/si14w/gnearline/ATAC_analysis/data/H7_DNase.distal.H3K27me3.bed  > ../txt/H7_DNase.distal.H3K27me3.txt
awk 'BEGIN{ print ("H3K9ac")} { print $NF}' /home/si14w/gnearline/ATAC_analysis/data/H7_DNase.distal.H3K9ac.bed  > ../txt/H7_DNase.distal.H3K9ac.txt

awk 'BEGIN{ print ("id")} { printf("input_%s\n",NR) }' ../data/H7_DNase.distal.bed > ../txt/H7_DNase_seqids.txt


paste ../txt/H7_DNase_seqids.txt ../txt/H7_DNase.distal.atac_score.txt ../txt/H7_DNase.distal.phylop_cons.txt ../txt/H7_DNase.distal.dist_to_gene.txt  ../txt/H7_DNase.distal.GC.txt ../txt/H7_DNase.distal.H3K4me1.txt ../txt/H7_DNase.distal.H3K4me3.txt ../txt/H7_DNase.distal.H3K9me3.txt ../txt/H7_DNase.distal.H3K27me3.txt ../txt/H7_DNase.distal.H3K9ac.txt ../txt/enhancer_distal_positive_negative_labels.txt |  awk 'BEGIN{print "id\tatac_score\tphylop_cons\tdist_to_gene\tgc\tca\tga\tgc_dinuc\tH3K4me1\tH3K4me3\tH3K9me3\tH3K27me3\tH3K9ac\tlabel"} { if(NR > 1) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"($6+$7)/$11"\t"($8+$9)/$11"\t"$10/$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$NF} }' > ../txt/H7_DNase_trainingset_withhismods.txt

awk '{ print $1":"$2"-"$3}' ../data/H7_DNase.distal.bed > ../txt/H7_DNase_distal.seqList
twoBitToFa ~/my/hg19.2bit ../txt/H7_DNase_peak_distal_sequences.fa -seqList=../txt/H7_DNase_distal.seqList


awk '{
                if (substr($1,1,1) != ">")
                {
                        printf("%s",$0)
                } else if (NR != 1) {
                        printf("\n")
                }
            } END { printf("\n") }' ../txt/H7_DNase_peak_distal_sequences.fa > ../txt/H7_DNase_peak_distal_sequence.only.txt


