module load cufflinks/2.2.1
echo """/home/si14w/gnearline/flu/stringtie_out/stringtie_out_mock12h.gtf
/home/si14w/gnearline/flu/stringtie_out/stringtie_out_6h.gtf
/home/si14w/gnearline/flu/stringtie_out/stringtie_out_12h.gtf
/home/si14w/gnearline/flu/stringtie_out/stringtie_out_24h.gtf
/home/si14w/gnearline/flu/stringtie_out/stringtie_out_48h.gtf""" > ../txt/cuffcompare_gtf_input_stringtie.txt
#cuffcompare -V -r ../txt/biocore_ucsc_plus_flu.gtf -i ../txt/cuffcompare_gtf_input_stringtie.txt -o ../txt/cuffcmp_stringtie
cuffmerge -p 12 -g ../txt/biocore_ucsc_plus_flu.gtf ../txt/cuffcompare_gtf_input_stringtie.txt -o ../txt/cuffmerge_stringtie
