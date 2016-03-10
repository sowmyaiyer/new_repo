module load star/2.4.2a
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /project/umw_flustore_schiffer/STAR_index_hg19_and_fluWithSnp_combined --genomeFastaFiles $HOME/gnearline/flu/bowtie_out/hg19_and_fluWithSnp_combined.fa --sjdbGTFfile $HOME/gnearline/flu/txt/biocore_ucsc_plus_flu.gtf
