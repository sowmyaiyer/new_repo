for 
samtools view ../bowtie_out/112213flucap_C.sorted.bam | grep "B59FULL" | awk '
        {
                if ($4 == 1 && index($6,"S") > 0)
                {
                        split($6,arr,"S");
                        if ( arr[1] ~ /^[0-9]*$/)
                        printf(">%s\n%s\n",,substr($10,1,arr[1]))
                }
        }' > ../txt/softClipLengths_C.txt
