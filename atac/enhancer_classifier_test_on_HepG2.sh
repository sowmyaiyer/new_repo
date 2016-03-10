#bedtools slop -i /project/umw_garberlab/tabakb/Maehr_AtacSeq/HumanES/22JUL14_PE50_C4MU1AC-B/Sample_HI_ATAC-Seq/Alignments/macs/HI_ATAC-Seq_NoIndex_L002.Merged.cutadapt.sorted_peaks.bed -b 50 -g ~/TR/hg19.genome >  H1_ATAC_macs2.bed
#bedtools slop -i /home/si14w/nearline/enhancer_predictions/gencodeV19_TSS.bed -l 1000 -r 500 -s -g ~/TR/hg19.genome > GencodeV19.promoter.bed


bedtools intersect -a /home/si14w/scratch/HepG2.DNase.peaks.narrowPeak -b  GencodeV19.promoter.bed -v > HepG2_DNase_macs2.distal.bed

awk '{ print $1"\t"$2"\t"$3"\t"NR}' ~/scratch/HepG2_DNase_macs2.distal.bed | bigWigAverageOverBed ~/farline/hot_mar2012/hg19.100way.phyloP100way.bw stdin ~/nearline/RECYCLE_BIN/out.tab -bedOut=/home/si14w/scratch/HepG2_DNase_macs2.distal.cons.bed

bedtools intersect -a HepG2_DNase_macs2.distal.bed -b HepG2.H3K27ac.broadPeak -wao  | bedtools groupby -i stdin -c 20 -o max | awk 'BEGIN{print "label"}{ if  ($NF == 0) label=-1; else label=1; print label}' > HepG2_test_labels.txt
bedtools nuc -fi ~/my/hg19.fa -bed HepG2_DNase_macs2.distal.bed | awk 'BEGIN{print ("seqLen")}{ if(NR > 1) print $19}' >  HepG2_DNase_macs2.distal.seqLen.txt
bedtools nuc -fi ~/my/hg19.fa -bed HepG2_DNase_macs2.distal.bed | awk 'BEGIN{print ("at""\t""gc")}{ if(NR > 1) print $11"\t"$12}' > HepG2_DNase_macs2.distal.AT_and_GC.txt
bedtools nuc -fi ~/my/hg19.fa -bed HepG2_DNase_macs2.distal.bed -pattern "CA" -C | awk 'BEGIN{ print ("ca" ) } { if (NR > 1) { printf("%.2f\n", $NF)} }' > HepG2_DNase_macs2.distal.CA.txt
bedtools nuc -fi ~/my/hg19.fa -bed HepG2_DNase_macs2.distal.bed -pattern "TG" -C | awk 'BEGIN{ print ("ca_rev" ) }{ if (NR > 1) { printf("%.2f\n", $NF)} }' > HepG2_DNase_macs2.distal.CA_rev.txt
bedtools nuc -fi ~/my/hg19.fa -bed HepG2_DNase_macs2.distal.bed -pattern "GA" -C | awk 'BEGIN{ print ("ga" ) } { if (NR > 1) { printf("%.2f\n", $NF)} }' > HepG2_DNase_macs2.distal.GA.txt
bedtools nuc -fi ~/my/hg19.fa -bed HepG2_DNase_macs2.distal.bed -pattern "TC" -C | awk 'BEGIN{ print ("ga_rev" ) } { if (NR > 1) { printf("%.2f\n", $NF)} }' > HepG2_DNase_macs2.distal.GA_rev.txt
bedtools nuc -fi ~/my/hg19.fa -bed HepG2_DNase_macs2.distal.bed -pattern "GC" -C | awk 'BEGIN{ print ("gc_dinuc" ) } { if (NR > 1) { printf("%.2f\n", $NF)} }' > HepG2_DNase_macs2.distal.GC.txt

awk 'BEGIN{ print ("atac_score" ) }{ print $7 }' HepG2_DNase_macs2.distal.bed > HepG2_DNase_macs2.distal.atac_score.txt
awk 'BEGIN{ print ("phylop_cons" ) }{ print $NF}' /home/si14w/scratch/HepG2_DNase_macs2.distal.cons.bed >  HepG2_DNase_macs2.distal.phylop_cons.txt
bedtools closest -a HepG2_DNase_macs2.distal.bed -b /home/si14w/nearline/enhancer_predictions/gencodeV19_TSS.bed -d -t "first" > HepG2_DNase_macs2_closest_tss.bed
awk 'BEGIN{ print ("dist_to_gene")} { print $NF}' HepG2_DNase_macs2_closest_tss.bed  > HepG2_DNase_macs2.distal.dist_to_gene.txt

awk 'BEGIN{ print ("id")} { printf("input_%s\n",NR) }' HepG2_DNase_macs2.distal.bed > HepG2_DNase_macs2_seqids.txt

paste HepG2_DNase_macs2_seqids.txt HepG2_DNase_macs2.distal.atac_score.txt HepG2_DNase_macs2.distal.phylop_cons.txt HepG2_DNase_macs2.distal.dist_to_gene.txt  HepG2_DNase_macs2.distal.AT_and_GC.txt HepG2_DNase_macs2.distal.CA.txt HepG2_DNase_macs2.distal.CA_rev.txt HepG2_DNase_macs2.distal.GA.txt HepG2_DNase_macs2.distal.GA_rev.txt HepG2_DNase_macs2.distal.GC.txt HepG2_DNase_macs2.distal.seqLen.txt enhancer_distal_positive_negative_labels.txt | awk 'BEGIN{print "id\tatac_score\tphylop_cons\tdist_to_gene\tat\tgc\tca\tga\tgc_dinuc\tlabel"} { if(NR > 1) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"($7+$8)/$(NF-1)"\t"($9+$10)/$(NF-1)"\t"$11/$(NF-1)"\t"$NF} }' > HepG2_DNase_macs2_trainingset.txt

awk '{ print $1":"$2"-"$3}' HepG2_DNase_macs2.distal.bed > HepG2_DNase_macs2_distal.seqList
twoBitToFa ~/my/hg19.2bit HepG2_DNase_macs2_peak_distal_sequences.fa -seqList=HepG2_DNase_macs2_distal.seqList


awk '{
                if (substr($1,1,1) != ">")
                {
                        printf("%s",$0)
                } else if (NR != 1) {
                        printf("\n")
                }
            } END { printf("\n") }' HepG2_DNase_macs2_peak_distal_sequences.fa > HepG2_DNase_macs2_peak_distal_sequence.only.txt

# Now get frequencies of kmer normalized by sequence length
##javac utils/KmerCounter.java
##java utils.KmerCounter $HOME/scratch/H1_ATAC_macs2_peak_distal_sequence.only.txt $HOME/scratch/H1_ATAC_macs2_distal_enhancer_positive_and_negative.counts.txt $HOME/scratch/kmers_6.txt # output - kmer counts

