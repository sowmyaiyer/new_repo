#!/bin/bash

module load samtools/0.0.19
files=$1 #bam file of flu alignment

echo "forward"
#forward reads by segments
samtools view -F 16 $files | awk ' ( $3 == "seg1_B59FULL" && $4 < 150 ){print $0}' | wc -l
samtools view -F 16 $files | awk ' ( $3 == "seg2_B59FULL" && $4 < 150 ){print $0}' | wc -l
samtools view -F 16 $files | awk ' ( $3 == "seg3_B59FULL" && $4 < 150 ){print $0}' | wc -l
samtools view -F 16 $files | awk ' ( $3 == "seg4_B59FULL" && $4 < 150 ){print $0}' | wc -l
samtools view -F 16 $files | awk ' ( $3 == "seg5_B59FULL" && $4 < 150 ){print $0}' | wc -l
samtools view -F 16 $files | awk ' ( $3 == "seg6_B59FULL" && $4 < 150 ){print $0}' | wc -l
samtools view -F 16 $files | awk ' ( $3 == "seg7_B59FULL" && $4 < 150 ){print $0}' | wc -l
samtools view -F 16 $files | awk ' ( $3 == "seg8_B59FULL" && $4 < 150 ){print $0}' | wc -l

echo "reverse"
#reverse reads by segments
samtools view -f 16 $files | awk ' ( $3 == "seg1_B59FULL" && $4 > 2295 ){print $0}' | wc -l
samtools view -f 16 $files | awk ' ( $3 == "seg2_B59FULL" && $4 > 2295 ){print $0}' | wc -l
samtools view -f 16 $files | awk ' ( $3 == "seg3_B59FULL" && $4 > 2190 ){print $0}' | wc -l
samtools view -f 16 $files | awk ' ( $3 == "seg4_B59FULL" && $4 > 1729 ){print $0}' | wc -l
samtools view -f 16 $files | awk ' ( $3 == "seg5_B59FULL" && $4 > 1519 ){print $0}' | wc -l
samtools view -f 16 $files | awk ' ( $3 == "seg6_B59FULL" && $4 > 1416 ){print $0}' | wc -l
samtools view -f 16 $files | awk ' ( $3 == "seg7_B59FULL" && $4 > 984 ){print $0}' | wc -l
samtools view -f 16 $files | awk ' ( $3 == "seg8_B59FULL" && $4 > 844 ){print $0}' | wc -l







