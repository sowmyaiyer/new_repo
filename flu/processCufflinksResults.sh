module load cufflinks/2.2.1
echo """/home/si14w/gnearline/flu/cufflinks_out_mock12h/transcripts.gtf
/home/si14w/gnearline/flu/cufflinks_out_6h/transcripts.gtf
/home/si14w/gnearline/flu/cufflinks_out_12h/transcripts.gtf
/home/si14w/gnearline/flu/cufflinks_out_24h/transcripts.gtf
/home/si14w/gnearline/flu/cufflinks_out_48h/transcripts.gtf""" > ../txt/cuffcompare_gtf_input.txt
cuffcompare -V -r ../txt/biocore_ucsc_plus_flu.gtf -i ../txt/cuffcompare_gtf_input.txt -o ../txt/cuffcmp
