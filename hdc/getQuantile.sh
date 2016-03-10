quartiles_jurkat=`awk '{ print $NF}' ../txt/gene_tss_expression_Jurkat.bed | Rscript getQuantile.R`
echo quartiles_jurkat
echo $quartiles_jurkat
quartiles_cd14_mono=`awk '{ print $NF}' ../txt/gene_tss_expression_roadmap.bed | Rscript getQuantile.R`
echo quartiles_cd14_mono
echo $quartiles_cd14_mono
quartiles_human_dc=`awk '{ print $NF}' ../txt/gene_tss_expression_time_0h.bed | Rscript getQuantile.R`
echo quartiles_human_dc
echo $quartiles_human_dc
quartiles_mouse_dc=`awk '{ print $NF}' ../txt/mouse_gene_tss_expression_time_mouse_0h.bed | Rscript getQuantile.R`
echo quartiles_mouse_dc
echo $quartiles_mouse_dc
